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
import subprocess
import time
from io_tools.file_utils import natural_sort_key, print_flush as print


def retrieve_chrom_size(g_liftoverSV):
    """
    Read a .chain file to memorize chromosome sizes of the target build.
    Stores in g_liftoverSV["size_chrom_target"][chrom]
    """
    # chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
    print(f"[{time.strftime('%H:%M:%S')}] Loading the size of the chromosomes of the target build")
    chain_file_name = os.path.basename(g_liftoverSV["chain"])
    print(f"           (parsing {chain_file_name})")
    g_liftoverSV["size_chrom_target"] = {}
    with open(g_liftoverSV["chain"], "r") as f:
        for line in f:
            if line.startswith("chain"):
                split_line = line.strip().split()
                # chrom_before_lift = split_line[2]
                # size_before_lift = int(split_line[3])
                # g_liftoverSV["size_before_lift"][chrom_before_lift] size_before_lift
                chrom_target = split_line[7]
                size_chrom_target = int(split_line[8])
                g_liftoverSV["size_chrom_target"][chrom_target] = size_chrom_target
  
    if g_liftoverSV["verbose"]:
        print("\n--verbose-- Size of the chromosomes of the target build:")
        n = 0
        n_up=len(g_liftoverSV["size_chrom_target"])-5
        for chrom, size in sorted(g_liftoverSV["size_chrom_target"].items(),
                                    key=lambda item: natural_sort_key(item[0])):
            n += 1
            if n <= 5 or n >= n_up: 
                print(f"{chrom}: {size}")
            if n == 6 and n_up > 6:
                    print(f"...")



def check_ref_fasta_seq(g_liftoverSV):
    """
    Check if chain file and reference FASTA file are coherent (with or without 'chr').
    """
    # chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2

    print(f"[{time.strftime('%H:%M:%S')}] Checking the ref_fasta_seq file")

    # Check the chain file
    chain_status = None
    with open(g_liftoverSV["chain"]) as f:
        for line in f:
            if line.startswith("chain"):
                chain_status = "with" if "chr" in line.split()[7] else "without"
                break

    # Check the ref_fasta_seq file
    fasta_status = None
    with open(g_liftoverSV["ref_fasta_seq"]) as f:
        for line in f:
            if line.startswith(">"):
                fasta_status = "with" if "chr" in line else "without"
                break

    if chain_status != fasta_status:
        print("\nIncoherence:")
        print("############")
        print(f"- chain: {g_liftoverSV['chain']}")
        print(f"     => Contig names of the target build: {chain_status} the 'chr' prefix")
        print(f"- ref_fasta_seq: {g_liftoverSV['ref_fasta_seq']}")
        print(f"     => Contig names of the target build: {fasta_status} the 'chr' prefix")
        print("\nPlease, check your ref_fasta_seq option")
        print("\nExit")
        sys.exit(2)


def is_multi_allelic(g_liftoverSV):
    """
    Check if the VCF input file contains multi-allelic lines.
    """

    print(f"[{time.strftime('%H:%M:%S')}] Ensuring that the input VCF contains only biallelic variants")

    input_file = g_liftoverSV["input_file"]

    if input_file.endswith(".gz"):
        cmd = f"zcat {input_file} | grep -v ^# | cut -f 4-5 | grep -c ,"
    else:
        cmd = f"grep -v ^# {input_file} | cut -f 4-5 | grep -c ,"

    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        output = result.stdout.strip()
        n_multiallelic_line = output[0] if output else "" 
        if n_multiallelic_line != "0":
            print("####################################################################################")
            print("Please split the multi-allelic lines of the VCF input file before to run liftoverSV")
            print("Exit without error.")
            print("####################################################################################")
            sys.exit(0)
    except subprocess.CalledProcessError:
        pass

