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

import gzip
import functools
import builtins
import gzip
import re
from itertools import islice



print_flush = functools.partial(builtins.print, flush=True)


def is_an_empty_vcf_file(vcf_file: str) -> bool:
    """
    Check if a VCF file (.vcf or .vcf.gz) is empty (ignoring header lines),
    or if it isn't a VCF file based on its extension.

    Returns:
        True if the file is empty or not a VCF, False otherwise.
    """
    # Quick check on the file extension
    if not re.search(r"\.vcf(\.gz)?$", vcf_file, re.IGNORECASE):
        return True

    # Determine whether to use gzip or normal open
    open_func = gzip.open if vcf_file.endswith(".gz") else open

    try:
        with open_func(vcf_file, "rt") as f:  # 'rt' = text mode
            for line in f:
                line = line.strip()
                # Ignore empty lines and header lines
                if line and not line.startswith("#"):
                    return False  # Found at least one data line
    except FileNotFoundError:
        # File doesn't exist → treat as empty
        return True
    except Exception as e:
        # Other read errors → treat as non-VCF/empty
        print(f"[WARNING] Could not read file {vcf_file}: {e}")
        return True

    # No data lines found → file is empty
    return True



def check_vcf_variant_line_format(vcf_path):
    """
    Check:
     - whether the #CHROM line exist
     - whether variant lines in a VCF file have the same number of fields
        as defined by the #CHROM header line.
     - wether empty variant line exists 

    Returns:
      "OK" if the VCF is well-formatted; otherwise, returns the error message.
    """
    
    chrom_header = None
    expected_field_number = None

    with open_any_text_file(vcf_path) as f:
        for line_number, line in enumerate(f, start=1):

            # Skip metadata
            if line.startswith("##"):
                continue

            # Extract #CHROM line
            if line.startswith("#CHROM"):
                chrom_header = line.rstrip("\n")
                expected_field_number = len(chrom_header.split("\t"))
                continue

            # Skip if #CHROM not found yet
            if chrom_header is None:
                return "Header line starting with #CHROM not found in the VCF!"

            # Stop if we reach an empty line
            if not line.strip():
                return "Empty line present in the VCF"

            # Process variant line
            fields = line.rstrip("\n").split("\t")
            num_fields = len(fields)

            if num_fields != expected_field_number:
                return f"Incorrect number of fields: {expected_field_number} expected, {num_fields} on line {line_number}"
                
    return "OK"




# =========================
# Check for "chr" prefix
# =========================

# Check for a VCF file (.vcf or .vcf.gz) or for a chain file
############################################################
# => Return "with" for a file with a "chr" prefix
# => Return "without" for a file without a "chr" prefix
#
# Header Lines of a chain file:
# chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
# score -- chain score
# tName -- chromosome (reference/target sequence)
# tSize -- chromosome size (reference/target sequence)
# tStrand -- strand (reference/target sequence)
# tStart -- alignment start position (reference/target sequence)
# tEnd -- alignment end position (reference/target sequence)
# qName -- chromosome (target sequence)
# qSize -- chromosome size (target sequence)
# qStrand -- strand (target sequence)
# qStart -- alignment start position (target sequence)
# qEnd -- alignment end position (target sequence)
# id -- chain ID
def file_with_chr(file_to_check):
    """
    Check if a VCF or chain file contains 'chr' prefixes.
    Returns 'with' or 'without'.
    """
    if file_to_check.endswith(".vcf.gz"):
        f = gzip.open(file_to_check, "rt")
        file_type = "vcf"
    elif file_to_check.endswith(".vcf"):
        f = open(file_to_check, "r")
        file_type = "vcf"
    elif file_to_check.endswith(".chain"):
        f = open(file_to_check, "r")
        file_type = "chain"
    else:
        return None

    with f:
        for line in f:
            line = line.strip()
            if file_type == "vcf":
                if line.startswith("#"):
                    continue
                return "with" if line.startswith("chr") else "without"
            elif file_type == "chain":
                if line.startswith("chain"):
                    chrom = line.split()[2]
                    return "with" if "chr" in chrom else "without"
    return None



# Natural sorting 
#################
# chromosomes = ["chr1", "chr10", "chr2", "chrX"]
# sorted_chromosomes = sorted(chromosomes, key=natural_sort_key)
# print(sorted_chromosomes)
# => Output: ['chr1', 'chr2', 'chr10', 'chrX']
def natural_sort_key(s):
    """Return a key for natural sorting (like Tcl -dictionary)."""
    # split into list of ints and non-ints
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]


def print_head(file_path: str, n: int = 5):
    """
    Print the first `n` lines of a file.
    
    Works for small or very large files without loading the entire file into memory.
    
    Args:
        file_path (str): Path to the file to read.
        n (int): Number of lines to print. Default is 5.
    """
    try:
        with open(file_path, "r") as f:
            for line in islice(f, n):
                print(line, end="")
    except FileNotFoundError:
        print(f"Error: file '{file_path}' not found.")


def open_any_text_file(path):
    """
    Open a text file transparently, whether it is plain text or gzip-compressed.

    This function automatically detects gzip compression by checking the magic
    bytes (0x1F, 0x8B), so it works regardless of the file extension. It is intended
    for reading text-based files (VCF, BED, TSV, CSV, config files, etc.) that may
    be distributed either uncompressed or compressed with gzip.

    Args:
        path (str): Path to the file to open.

    Returns:
        file object: A text-mode file handle. 
        If the file is gzip-compressed, a gzip.open() handle is returned; otherwise, a standard open().

    Raises:
        OSError: If the file cannot be opened or read.
    """
    # Detect gzip by reading the first two bytes of the file
    with open(path, "rb") as test:
        start = test.read(2)

    # Magic bytes for gzip compression = 1F 8B
    if start == b"\x1f\x8b":
        return gzip.open(path, "rt")   # gzip-compressed text file
    else:
        return open(path, "rt")        # plain-text file



