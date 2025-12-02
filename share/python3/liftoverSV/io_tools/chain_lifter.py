
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

from typing import List, Tuple, Optional
from pyliftover import LiftOver


class ChainLifter:
    """
    A utility class for performing genomic coordinate liftOver operations
    using UCSC .chain files.

    This version ignores the strand information and only returns the
    new vcf_chromosome and 1-based position for each input coordinate.

    Conversion:
        Input coordinates must be 1-based (VCF convention).
        Output coordinates are also 1-based.

    Example:
        chain = ChainLifter("hg19ToHg38.over.chain.gz")
        new_chrom, new_pos = chain.lift("chr1", 12345)
        results = chain.lift_many([("chr1", 12345), ("chr2", 67890)])
    """
    def __init__(self, chain_file):
        """Load the .chain file once (using pyliftover) and returns the LiftOver object."""
        self.lo = LiftOver(chain_file)

    def lift(self, vcf_chrom: str, vcf_pos: int) -> Tuple[Optional[str], Optional[int]]:
        """
        Lift a single genomic coordinate (vcf_chromosome, vcf_position).
        
        Args:
            vcf_chrom (str): Chromosome name in the source genome
            vcf_pos (int): 1-based coordinate in the source genome

        Returns:
            Tuple[str|None, int|None]: New chromosome and position
                                        or (None, None) if lift fails
        """
        # Convert 1-based VCF coordinate to 0-based for LiftOver
        bed_pos=vcf_pos-1

        # "convert_coordinate" method: 
        # - Input coordinates must be 0-based
        # - output coordinates are also 0-based.
        lifted = self.lo.convert_coordinate(vcf_chrom, bed_pos)
        if lifted:
            new_vcf_chrom, new_bed_pos, strand, _chain = lifted[0]
            # Convert back to 1-based
            new_vcf_pos=new_bed_pos+1
            return new_vcf_chrom, new_vcf_pos
        else:
            return None, None

    def lift_many(self, positions: List[Tuple[str, int]]) -> List[Tuple[Optional[str], Optional[int]]]:
        """
        Lift multiple genomic coordinates in batch.
        Args:
            positions (list[tuple[str, int]]): List of (vcf_chrom, pos) tuples.
        Returns:
            list[tuple[str|None, int|None]]:
        """
        return [self.lift(vcf_chrom, vcf_pos) for vcf_chrom, vcf_pos in positions]
