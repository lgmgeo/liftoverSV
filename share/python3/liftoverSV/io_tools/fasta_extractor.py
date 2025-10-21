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
from typing import List, Tuple
from pyfaidx import Fasta


class FastaExtractor:
    """
    Efficiently extract sequences from a reference FASTA file using pyfaidx.
    Automatically uses the .fai index if it exists, otherwise creates it
    and informs the user.
    
    Usage:
        extractor = FastaExtractor("/path/to/ref.fasta")
        seq = extractor.get_sequence("chr1", 1000, 1010)
    """

    def __init__(self, fasta_path: str):
        """
        Initialize the extractor and use pyfaidx for indexed access.

        Args:
            fasta_path (str): Path to the reference FASTA file
        """
        self.fasta_path = fasta_path
        fai_path = fasta_path + ".fai"
        index_was_missing = not os.path.exists(fai_path)

        # Check if .fai exists
        if index_was_missing:
            print(fai_path)
            print(f"[INFO] Index file '{fai_path}' does not exist. Creating it...")

        # Load FASTA with pyfaidx; creates index if needed
        self.extractor = Fasta(fasta_path)

        # Confirm index is ready
        if index_was_missing and os.path.exists(fai_path):
            print(f"[INFO] Index '{fai_path}' has been created.")
            

    def get_sequence(self, vcf_chrom: str, vcf_start: int, vcf_end: int) -> str:
        """
        Extract a sequence from the reference FASTA.

        Args:
            vcf_chrom (str): Chromosome/contig name
            vcf_start (int): 1-based  start position 
            vcf_end (int): 1-based end position (inclusive)

        Returns:
            str: DNA sequence from start to end
        
        Raises:
            ValueError: if chromosome not found in FASTA

        """
        # Convert to 0-based indexing for pyfaidx
        bed_start = vcf_start - 1

        try:
            return self.extractor[vcf_chrom][bed_start:vcf_end].seq
        except KeyError:
            raise ValueError(f"Chromosome '{vcf_chrom}' not found in FASTA.")


    def get_sequences_from_bed(self, bed_lines: List[Tuple[str, int, int]]) -> List[str]:
        """
        Extract multiple sequences from a list of BED-like coordinates.

        Args:
            bed_lines (list of tuples): Each tuple is (chrom, start, end)

        Returns:
            list of str: Extracted DNA sequences
        """
        sequences = []
        for chrom, start, end in bed_lines:
            sequences.append(self.get_sequence(chrom, start, end))
        return sequences
    
    
    
    