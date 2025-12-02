"""
liftoverSV 0.3.0_beta
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
import time
import tempfile
from io_tools.file_utils import open_any_text_file, print_flush as print
from io_tools.batch_writer import BatchWriter
from io_tools.chain_lifter import ChainLifter
from io_tools.fasta_extractor import FastaExtractor
from core.liftover_engine import LiftoverEngine, Variant
from core.header_tools import extract_header_ids
from workflow.output_writer import add_new_header_lines, sort_and_compress_the_lifted_vcf
from multiprocessing import Pool, cpu_count

def process_chunk(chunk, g_liftoverSV):
    """
    Process a single chunk of VCF lines by performing the liftover operation.

    This function is executed in parallel by worker processes.  
    For each chunk, it:
      - Initializes a local ChainLifter, FastaExtractor, and LiftoverEngine instance
      - Converts each VCF line into a Variant object
      - Applies the liftover transformation to each variant
      - Collects both lifted and unmapped variants
      - Returns the results along with metadata collected by the LiftoverEngine

    Parameters
    ----------
    chunk : list of tuple(int, str)
        A list of (line_number, vcf_line) pairs 
    g_liftoverSV : dict
   
    Returns
    -------
    tuple
        A tuple containing:
        - results: list of (lifted_variant, reason) pairs
        - S_SVlines_INFO  : set of INFO IDs extracted from lifted variants
        - S_SVlines_FORMAT: set of FORMAT IDs extracted from lifted variants
        - S_SVlines_FILTER: set of FILTER IDs extracted from lifted variants
        - S_lifted_contigs: set of contigs observed during liftover
        - case_counts     : dict with counts for different SV categories
        - n_mapped        : number of successfully lifted variants
        - n_unmapped      : number of unmapped variants
    """    
    # Load the chain file
    chain = ChainLifter(g_liftoverSV['chain'])

    # Load FASTA with pyfaidx: uses .fai index if available
    extractor = FastaExtractor(g_liftoverSV["ref_fasta_seq"])

    # Initialise the LiftoverEngine
    engine = LiftoverEngine(chain, extractor, g_liftoverSV["percent"])

    results = []
    for line_number, line in chunk:
        # Create a Variant object from the line
        variant = Variant.from_vcf_line(line, line_number)
        # Lift the variant
        lifted, reason = engine.lift_variant(variant, g_liftoverSV)
        # Memorize the results
        results.append((lifted, reason))

    return results, engine.S_SVlines_INFO, engine.S_SVlines_FORMAT, engine.S_SVlines_FILTER, engine.S_lifted_contigs, engine.case_counts, engine.n_mapped, engine.n_unmapped



# Source:
#########
# https://samtools.github.io/hts-specs/VCFv4.4.pdf
#
# END:
######
# The END of each allele is defined as:
#   - Non-symbolic alleles: POS + length of REF allele − 1.
#   - <INS> symbolic structural variant alleles: POS + length of REF allele − 1.
#   - <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles: POS + SVLEN.
# <*> symbolic allele: the last reference call position.
# END must be present for all records containing the <*> symbolic allele
#
# SVLEN:
########
# SVLEN is defined for INS, DUP, INV , and DEL symbolic alleles as the number of the inserted, duplicated, inverted, and deleted bases respectively.
# SVLEN is defined for CNV symbolic alleles as the length of the segment over which the copy number variant is defined.
# The missing value . should be used for all other ALT alleles, including ALT alleles using breakend notation.
# ==> For INS, the SVLEN is defined as the number of the inserted bases. No liftover needed!!

# INFO/CIPOS and INFO/CIEND:
############################
# Check and modify if needed the CIPOS/CIEND values in order to have:
#	=> POS-CIPOS > 0
#	=> END+CIEND < chrom_length
def write_the_lifted_vcf(g_liftoverSV):
    """
    Perform liftover of a structural variant VCF using a chain file.
    - Remove all header lines (will be updated later in "def add_new_header_lines")
    - Generates unique IDs for variants missing one ("lifted_from_l_<line_number>").
    - Writes mapped variants to output_file (gzipped) and
    - unmapped/dropped variants to unmapped_file (tab-delimited).
    
    Rules to write the lifted VCF:
    => Lift #CHROM, POS, REF, ALT, INFO/END and INFO/SVEND
    => Drop the SV if:
        Case 1: One or more required positions fail to lift
        - "Start" coordinate
        - "End" coordinate
        - "Last REF base" coordinate
        - "Square-bracketed ALT" coordinate
        - "Last base in REF before a deleted segment" coordinate

        Case 2: Lifted positions map to different chromosomes (except for translocations):
        - "Start" and "end" coordinates
        - "Start" and "last REF base" coordinates
        - "Start" and "square-bracketed ALT" coordinates
        - "Start" and "last base in REF before a deleted segment" coordinates

        Case 3: Reversed order between lifted positions:
        - "Start" and "end" coordinates
        - "Start" and "square-bracketed ALT" coordinates

        Case 4: Significant change in distance after liftover (default: >5% of SVLEN):
        - Distance between "start" and "end" coordinates changed significantly
        - Distance between "start" and "square-bracketed ALT" coordinates changed significantly

        Case 5: Complex or inconsistent REF/ALT sequences
        - REF contains '.' or '*' inside the sequence
        - ALT contains '.' or '*' inside the sequence
        - Deletion: ALT sequence not at the beginning of REF sequence (e.g., REF="ATTCTTG", ALT="TC")
        - Insertion: REF sequence not at the beginning of ALT sequence (e.g., REF="TC", ALT="ATTCTTG")
        - Insertion with single breakend: REF sequence not opposite to '.' in ALT (e.g., REF="G", ALT=".TTTTTTC")
        - Square-bracketed ALT badly formatted
        - The REF sequence differs from the original after liftover (see REF and ALT)                                                                                                                                                                    
    """

    # Create a temporary VCF output file for mapped variants unsort (no header lines)
    print(f"[{time.strftime('%H:%M:%S')}] Initializing a temporary VCF output file for unsorted mapped variants")
    tmp_output_file = tempfile.NamedTemporaryFile(delete=False, dir=g_liftoverSV["tmp_dir"], suffix=f".liftoverSV.tmp.vcf").name
    print(f"           {tmp_output_file}")
    tmp_out_writer = BatchWriter(tmp_output_file, g_liftoverSV)    

    # Initialize a BatchWriter for the unmapped file
    print(f"[{time.strftime('%H:%M:%S')}] Initializing the output unmapped file")
    output_file = g_liftoverSV['output_file']
    unmapped_file = re.sub(r"\.sort\.vcf.gz$", ".unmapped", output_file)
    if os.path.exists(unmapped_file):
        os.remove(unmapped_file)
    unmapped_writer = BatchWriter(unmapped_file, g_liftoverSV)
    print(f"           {unmapped_file}")

    # Input VCF (gzipped or not)
    input_file = g_liftoverSV['input_file']
    print(f"[{time.strftime('%H:%M:%S')}] Reading input VCF: {input_file}")

    # - Increment "vcf_line_number"
    # - Lift over SV
    ####################################################################
    # All "lifted" variants should have a unique identifier:
    # - Missing IDs (ID=".") are replaced with "lifted_from_<vcf_line_number>""
    # - Existing IDs are preserved.
    # ==> use the "vcf_line_number" variable (to later create "lifted_from_l_<line_number>" indices)
    print(f"[{time.strftime('%H:%M:%S')}] Lift over SV:")
    print(f"           Writing to {tmp_output_file}")
    print(f"           Writing to {unmapped_file}")
    at_least_1_SV_lifted = 0

    # Collect lines to parallelize
    L_chunks = []
    chunk = []

    # Initialize variables to update the header lines
    S_header_INFO, S_header_FORMAT, S_header_FILTER = set(), set(), set()
    S_SVlines_INFO, S_SVlines_FORMAT, S_SVlines_FILTER = set(), set(), set()
    S_lifted_contigs = set()
    # Counters for unmapped statistics
    n_mapped = 0
    n_unmapped = 0
    case_counts = {f"case{i}": 0 for i in range(1, 6)}
    
    with open_any_text_file(input_file) as f:
        for vcf_line_number, line in enumerate(f, 1):

            if line.startswith("#"):
                # Updade S_header_INFO, S_header_FORMAT and S_header_FILTER
                S_header_INFO, S_header_FORMAT, S_header_FILTER = extract_header_ids(line, S_header_INFO, S_header_FORMAT, S_header_FILTER)
            else:
                chunk.append((vcf_line_number, line))
                if len(chunk) >= g_liftoverSV["chunk_size"]:
                    L_chunks.append(chunk)
                    chunk = []                  
        if chunk:
            L_chunks.append(chunk)

    total_chunks = len(L_chunks)
    print(f"[{time.strftime('%H:%M:%S')}] Processing {total_chunks} chunks (target chunk size: {g_liftoverSV['chunk_size']} lines)")


    # Determine the number of workers to use: 
    # - value from g_liftoverSV if provided
    # - otherwise default to the number of CPU cores
    # WARNING: Ensure that the number of workers does not exceed the number of CPU cores
    #          to avoid unnecessary context switching and overhead.
    n_workers = g_liftoverSV.get("n_workers", cpu_count())
    if n_workers > cpu_count():
       n_workers = cpu_count() 
    
    # Prepare the tuples
    args = [(chunk, g_liftoverSV) for chunk in L_chunks]
    
    i_chunk = 1
    with Pool(n_workers) as pool:
        for (
            result_batch,
            Si_SVlines_INFO,
            Si_SVlines_FORMAT,
            Si_SVlines_FILTER,
            Si_lifted_contigs,
            i_case_counts,
            ni_mapped,
            ni_unmapped
        ) in pool.starmap(process_chunk, args):
            if g_liftoverSV["verbose"]:
                print(f"[{time.strftime('%H:%M:%S')}] Chunk {i_chunk}/{total_chunks}")
            i_chunk += 1

            for lifted_variant, reason in result_batch:
                if lifted_variant:
                    # Liftover successfull
                    tmp_out_writer.write(lifted_variant)
                    at_least_1_SV_lifted = 1 
                else:
                    # Unmapped
                    unmapped_writer.write(reason)

            S_SVlines_INFO.update(Si_SVlines_INFO)
            S_SVlines_FORMAT.update(Si_SVlines_FORMAT)
            S_SVlines_FILTER.update(Si_SVlines_FILTER)
            S_lifted_contigs.update(Si_lifted_contigs)
            n_mapped += ni_mapped
            n_unmapped += ni_unmapped
            case_counts = {k: case_counts.get(k, 0) + i_case_counts.get(k, 0) for k in set(case_counts) | set(i_case_counts)}

    
    # Close and flush the remaining lines:
    tmp_out_writer.close()
    unmapped_writer.close()

    print(f"[{time.strftime('%H:%M:%S')}] Liftover summary:")
    print(f"           * {n_mapped} mapped SV")
    print(f"           * {n_unmapped} unmapped SV")
    if n_unmapped:
        print(f"           (see {unmapped_file} for details)")

    # Map cases to their messages
    case_messages = {
        "case1": "SVs where one or more required positions failed to lift",
        "case2": "SVs where two positions (start, end, etc.) mapped to different chromosomes (except for translocations)",
        "case3": "SVs where lifted positions are in reverse order",
        "case4": "SVs where the distance between positions changed significantly after liftover (default: >5% of SVLEN)",
        "case5": "SVs with complex or inconsistent REF/ALT sequences"
    }        

    # Only print cases that are non-zero / truthy
    for case, message in case_messages.items():
        if case_counts.get(case, 0):
            print(f"             - {case_counts[case]} {message}")

    # Exit if no SV lifted
    if not at_least_1_SV_lifted:
        print(f"[{time.strftime('%H:%M:%S')}] Liftover completed successfully.")
        sys.exit(0)

    # Keep only new header lines (only for INFO, FORMAT, or FILTER) that are not already in the existing headers
    S_new_INFO   = list(S_SVlines_INFO - S_header_INFO)
    S_new_FORMAT = list(S_SVlines_FORMAT - S_header_FORMAT)
    S_new_FILTER = list(S_SVlines_FILTER - S_header_FILTER)

    # Add them to the VCF header
    # assembly / contig / FILTER / INFO / FORMAT
    # => creation of g_liftoverSV['output_file']
    add_new_header_lines(S_new_INFO, S_new_FORMAT, S_new_FILTER, S_lifted_contigs, input_file, g_liftoverSV)
   
    # Sort the output file
	######################
    sort_and_compress_the_lifted_vcf(tmp_output_file, g_liftoverSV)

