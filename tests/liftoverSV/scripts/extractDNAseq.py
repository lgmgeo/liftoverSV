#!/usr/bin/env python3

# Command line example:
# python3 ./extractDNAseq.py chr22 16848507 16848507 $ref_fasta_seq
# => g
# python3 ./extractDNAseq.py chr22 16848507 16848517 $ref_fasta_seq
# => gcatatccttg

import sys
import subprocess

def extract_dna_seq(vcf_chrom, vcf_start, vcf_end, ref_fasta):
    """
    Extract DNA sequence from a reference FASTA using 1-based VCF coordinates.

    Args:
        vcf_chrom (str): Chromosome (e.g., 'chr22')
        vcf_start (int): Start position (1-based)
        vcf_end (int): End position (1-based, inclusive)
        ref_fasta (str): Path to reference FASTA file

    Returns:
        str: DNA sequence, or empty string in case of error
    """
    # Convert 1-based VCF to 0-based BED start
    bed_start = vcf_start - 1

    # Prepare bedtools command
    # echo -e "chr start end" | bedtools getfasta -fi ref.fa -bed -
    bed_line = f"{vcf_chrom}\t{bed_start}\t{vcf_end}"
    cmd = ["bedtools", "getfasta", "-fi", ref_fasta, "-bed", "-"]

    try:
        # Run bedtools with BED coordinates on stdin
        result = subprocess.run(
            cmd,
            input=bed_line,
            text=True,
            capture_output=True,
            check=True
        )
        # Remove FASTA header lines (starting with ">") and empty lines
        seq_lines = [line.strip() for line in result.stdout.splitlines()
                     if line.strip() and not line.startswith(">")]
        return "".join(seq_lines)
    except subprocess.CalledProcessError as e:
        # Print command and error message for debugging
        print("Error running command:", " ".join(cmd), file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        return ""

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} <chrom> <start> <end> <ref_fasta>", file=sys.stderr)
        sys.exit(1)

    chrom = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    ref_fasta = sys.argv[4]

    sequence = extract_dna_seq(chrom, start, end, ref_fasta)
    print(sequence)


