"""
liftoverSV 0.3.1_beta
=====================

Copyright (C) 2024-current Veronique Geoffroy (veronique.geoffroy@inserm.fr)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import re
import subprocess
import shutil
import argparse
import tempfile
from io_tools.file_utils import is_an_empty_vcf_file, file_with_chr, check_vcf_variant_line_format, print_flush as print
from functools import partial
from core.constants import CHUNK_SIZE, N_WORKERS


def valid_vcf_input_file(vcf_input_file):
    """
    Check if the VCF input file:
    - exists
    - has .vcf or .vcf.gz extension
    - contains at least 1 SV line (not empty)
    - contains the #CHROM line
    - has the same number of fields in the variant lines as defined by the #CHROM header line.
    - empty variant line exists 
    """
    # Check if the file exists
    if not os.path.isfile(vcf_input_file):
        # Normal error -> argparse will handle and exit with code 2
        raise argparse.ArgumentTypeError(f"File '{vcf_input_file}' does not exist")

    # Check if the file has a valid VCF extension
    if not (vcf_input_file.endswith(".vcf") or vcf_input_file.endswith(".vcf.gz")):
        raise argparse.ArgumentTypeError(f"File '{vcf_input_file}' must have .vcf or .vcf.gz extension")

    # Special case: empty VCF file -> print message and exit with code 0
    if is_an_empty_vcf_file(vcf_input_file):
        print(f"File '{vcf_input_file}' is empty, no SV to lift")
        sys.exit(0)

    # Check if the input VCF file contains multi-allelic lines
    message = check_vcf_variant_line_format(vcf_input_file)
    if message != "OK":
        print(message)
        sys.exit(1)

    # If all checks pass, return the vcf_input_file
    return vcf_input_file
    

def valid_chain_file(chain_file):
    """
    Validate the liftover chain file:
    - must exist
    - must not be empty
    - must end with '.chain'
    - must be consistent with the input VCF file (chr prefix)
    """

    # Check if file exists
    if not os.path.exists(chain_file):
        print("\n############################################################################")
        print(f"Bad value for the --chain option, file does not exist ({chain_file})")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # Check if the provided "chain_file" value is empty
    if not chain_file:
        print("\n############################################################################")
        print("liftoverSV needs in argument the path of your chain file (--chain ...)")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # Check extension
    if not re.search(r"\.chain$", chain_file, re.IGNORECASE):
        print("\n############################################################################")
        print(f"Bad option value: --chain = {chain_file}")
        print("Extension file should be '.chain'")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)
    
    # If all checks pass, return the chain_file
    return chain_file


def valid_ref_fasta(ref_fasta):
    """
    Validate the reference FASTA file:
    - file must not be empty
    - file must have .fasta or .fa extension
    - file must exist
    """
    # Check if a path was provided
    if not ref_fasta:
        print("\n############################################################################")
        print("liftoverSV needs in argument the path of your ref-fasta-seq file (--ref-fasta-seq ...)")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # Check valid extension
    if not re.search(r"\.(fasta|fa)$", ref_fasta, re.IGNORECASE):
        print("\n############################################################################")
        print(f"Bad option value: --ref-fasta-seq = {ref_fasta}")
        print("Extension file should be '.fasta' or '.fa'")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # Check if the file exists
    if not os.path.exists(ref_fasta):
        print("\n############################################################################")
        print(f"Bad value for the --ref-fasta-seq option, file does not exist ({ref_fasta})")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # If all checks pass, return the ref_fasta
    return ref_fasta


def valid_tool_path(tool_path, tool_name):
    """
    Validate that the given tool is installed and executable.
    - If a full path is provided, it must exist
    - If only a command name is provided, it must be found in $PATH
    - Validate that a CLI tool exists and prints 'usage' or 'help' when run without arguments
    """
        
    # Resolve "full path" / "command name" in PATH
    resolved_path = shutil.which(tool_path)
    if resolved_path is None:
        print(f"\nError: {tool_name} not found in PATH ('{tool_path}').")
        sys.exit(2)

    # Try running the tool
    try:
        # Run the tool without arguments, capture stdout and stderr
        result = subprocess.run(
            [resolved_path],
            stdout=subprocess.PIPE,  # capture stdout
            stderr=subprocess.STDOUT,  # redirect stderr to stdout
            text=True,                # return string instead of bytes
            timeout=5                  # optional: avoid hanging
        )
 
        # Check if 'usage' or 'help' appears in output
        output = result.stdout.lower()
        if "usage" not in output and "help" not in output:
            print(f"\nError: {tool_name} does not seem valid ('{tool_path}').")
            sys.exit(2)

    except Exception:
        print(f"\nError: Cannot execute {tool_name} ('{tool_path}'). {str(e)}")
        sys.exit(2)

    return resolved_path

       
def valid_percent(value):
    """
    Validate the percent argument:
    - must be a float
    - must be in the [0, 1] range
    """
    try:
        percent = float(value)
    except ValueError:
        print(f"\nError: --percent must be a float, got '{value}'")
        sys.exit(2)

    if not (0 <= percent <= 1):
        print("\n############################################################################")
        print(f"Bad option value: --percent = {percent}")
        print("Should be in the [0-1] range values, default = 0.05")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    return percent


def valid_output_base_name(output_base_name):
    """
    Validate and normalize the output-file path.
    
    Checks:
    - File name is provided
    - If needed, remove ".sort.vcf" or ".vcf" and remove ".gz"
    - Output directory exists if a path is provided (else, if not provided, the output path is set to ".")

    Returns the normalized output_base_name
    """
    if not output_base_name:
        print("\n############################################################################")
        print("liftoverSV needs in argument the path of your output-file file (--output-file ...)")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

    # Normalize the output file name
    normalized_base_name = re.sub(r"\.gz$", "", output_base_name)
    normalized_base_name = re.sub(r"\.sort\.vcf$", "", normalized_base_name)
    normalized_base_name = re.sub(r"\.vcf$", "", normalized_base_name)

    # Return normalized file path
    return normalized_base_name

    
def configure_liftover_sv(argv, g_liftoverSV):
    """
    Configure liftoverSV options from argv.
    """

    # Creation of the parser
    ########################
    parser = argparse.ArgumentParser(
        description="Lift over structural variations (SV) between genome builds",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Definition of the arguments
    #############################

    # ───────────────────────────────────────────
    # 0) HELP & VERSION
    # (argparse automatically adds -h/--help)
    # ───────────────────────────────────────────
    parser.add_argument(
        "-V", "--version",
        action="version",
        version=f"liftoverSV {g_liftoverSV['version']}",
        help="show program's version number and exit"
    )    

    # ───────────────────────────────────────────
    # 1) INPUT FILES
    # ───────────────────────────────────────────
    group_input = parser.add_argument_group("Input files")

    group_input.add_argument(
        "-c", "--chain",
        type=valid_chain_file,
        required=True,
        metavar="<File>",
        help="""the liftover chain file
see https://genome.ucsc.edu/goldenPath/help/chain.html for a description of chain files
see http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files
required"""
    )

    group_input.add_argument(
        "-i", "--input-file", dest="input_file",
        metavar="<File>",
        required=True,
        type=valid_vcf_input_file,
        help="""the SV VCF input file
gzipped VCF file is supported
multi-allelic lines are not allowed
required"""
    )

    group_input.add_argument(
         "-r", "--ref-fasta-seq", dest="ref_fasta_seq",
         type=valid_ref_fasta,
         required=True,
         metavar="<File>",
         help="""the reference sequence (fasta) for the TARGET genome build (i.e. the new one after the liftover)
required"""         
    )

    # ───────────────────────────────────────────
    # 2) OUTPUT OPTIONS
    # ───────────────────────────────────────────
    group_output = parser.add_argument_group("Output options")

    group_output.add_argument(
        "-d", "--output-dir", dest="output_dir",
        type=str, 
        metavar="<Dir>",
        help="""the liftover SV VCF output directory
default: current directory"""
    )

    group_output.add_argument(
        "-o", "--output-base-name", dest="output_base_name",
        type=partial(valid_output_base_name),
        required=True,
        metavar="<File>",
        help="""Base name for output (generates FILE.sort.vcf.gz and FILE.unmapped)
required"""
    )

    # ───────────────────────────────────────────
    # 3) PERFORMANCE / PARALLELIZATION
    # ───────────────────────────────────────────
    group_perf = parser.add_argument_group("Performance")

    group_perf.add_argument(
    "-w", "--n-workers", type=int, default=N_WORKERS,
    help="""number of parallel worker processes to use for liftover.
increasing the number of workers can speed up processing on multi-core machines.
default: 8"""
    )

    group_perf.add_argument(
        "-z", "--chunk-size", type=int, default=CHUNK_SIZE,
        help="""number of VCF lines to process per chunk.
processing the VCF in chunks reduces memory usage and enables parallel liftover.
default: 50000"""
    )

    # ───────────────────────────────────────────
    # 4) BEHAVIORAL PARAMETERS
    # ───────────────────────────────────────────
    group_behavior = parser.add_argument_group("Behavior")

    group_behavior.add_argument(
         "-p", "--percent", 
         type=valid_percent,
         default=0.05,
         metavar="<float>",
         # Note: use '%%' instead of '%' in help strings.
         # Argparse interprets '%' as a format specifier (e.g., %(default)s),
         # so a single '%' in the text will raise a ValueError during help display.
         help="""variation in length authorized for a lifted SV (e.g. difference max between both SVLENs < 5%%)
default value: 0.05""" 
    )

    group_behavior.add_argument(
        "-T", "--tmp-dir", dest="tmp_dir",
        type=str,
        metavar="<Dir>",
        help="""Directory where temporary files will be created.
    If not provided, the system default temporary directory is used."""
    )

    group_behavior.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="enable verbose output"
    )

    # Parsing of the arguments
    ##########################
    # WARNING:
    # Argparse converts option names to valid Python identifiers:
    # - long names become lowercase
    # - hyphens '-' are replaced with underscores '_'
    # Access the value via args.option_name, e.g., input-file >> args.input_file
    args = parser.parse_args()

    # Completion of the g_liftoverSV dictionary
    ###########################################
    g_liftoverSV.update(vars(args))

    # Check tmp_dir
    ###############
    # Determine tmp_dir
    if g_liftoverSV["tmp_dir"] is None:
        g_liftoverSV["tmp_dir"] = tempfile.gettempdir()  # default system tmp
    else:
        # Ensure directory exists
        if not os.path.isdir(g_liftoverSV["tmp_dir"]):
            raise ValueError(f"Temporary directory does not exist: {g_liftoverSV['tmp_dir']}")
    
    # Determine output_dir if not given in argument
    ###############################################
    if g_liftoverSV["output_dir"] == None:
        if "/" in g_liftoverSV["output_base_name"]:
            output_dir = os.path.dirname(g_liftoverSV["output_base_name"])
            if not os.path.exists(output_dir):
                output_dir = "."    
        else:
            output_dir = "."
    # Store output_dir in global dictionary
    g_liftoverSV["output_dir"] = output_dir

    # Determine output_file
    #######################
    g_liftoverSV["output_file"] = g_liftoverSV["output_base_name"] + ".sort.vcf.gz"

    # Cross-check "chr" prefix consistency between input-file and chain
    ##################################################################
    #chain_file = g_liftoverSV.get("chain", "")
    if file_with_chr(args.input_file) != file_with_chr(args.chain):
        print("\n############################################################################")
        print("Bad chain file:")
        print(f"- input file {file_with_chr(args.input_file)} prefix 'chr' ({args.input_file})")
        print(f"- chain file {file_with_chr(args.chain)} prefix 'chr' ({args.chain})")
        print("Exit with error.")
        print("############################################################################\n")
        sys.exit(2)

