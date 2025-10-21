"""
liftoverSV 0.2.0_beta
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

"""
variant.py

Defines the Variant class which represents a single variant record in a VCF file.

Responsibilities:
- Parse a VCF line into a structured object
- Automatically handle missing IDs
- Provide a method to reconstruct a VCF line after lifting or modifications
"""

from dataclasses import dataclass
from typing import List, Optional
from core.svtype_normalization import normalize_svtype
from io_tools.file_utils import print_flush as print



@dataclass
class Variant:
    """
    Variant represents a single variant (SV or SNV) from a VCF file.

    Attributes:
        chrom (str): #CHROM
        pos (int): POS (1-based)
        id (str): ID (auto-generated if missing)
        ref (str): REF
        alt (str): ALT
        qual (str): QUAL
        filter (str): FILTER
        info (str): INFO 
        format (str): FORMAT (optional)
        samples (List[str]): Genotype/sample fields (optional)
        line_number (int): Line number in the input VCF (used for auto-generated ID)
    """
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: str
    filter: str
    info: str
    var_head: str  
    svtype: Optional[str]
    format: str = ""
    ref_length: Optional[int] = None
    last_nt_coord: Optional[int] = None   
    samples: Optional[List[str]] = None

    line_number: Optional[int] = None


    # Lifted fields (initially unknown)
    lifted_chrom: Optional[str] = None
    lifted_pos: Optional[int] = None
    lifted_id: Optional[str] = None
    lifted_ref: Optional[str] = None
    lifted_last_nt_chrom: Optional[int] = None
    lifted_last_nt_coord: Optional[int] = None
    lifted_alt: Optional[str] = None
    lifted_info: Optional[str] = None

    @classmethod
    def from_vcf_line(cls, line: str, line_number: int):
        """
        Parse a VCF line and create a Variant object.

        Args:
            line (str): A line from a VCF file (tab-separated)
            line_number (int): Line number in the VCF file
                               (used for generating a unique ID if missing)

        Returns:
            Variant: An instance of the Variant class

        Notes:
            - If the ID field in VCF is ".", it is replaced by "lifted_from_l_<line_number>"
            - Handles optional FORMAT and sample columns
        """
        fields = line.rstrip().split("\t")  # Remove trailing newline and split by tab

        chrom = fields[0]  
        pos = int(fields[1])  # 1-based position

        # All lifted variants have a unique identifier (ID):
        # - Missing IDs (ID=".") are replaced with "lifted_from_<line_from_input_file>"
        # - Existing IDs are preserved.
        raw_id = fields[2]  
        if raw_id == ".":
            # Auto-generate a unique ID for variants with missing ID
            id = f"lifted_from_l_{line_number}"
        else:
            id = raw_id

        ref = fields[3]  # REF allele
        alt = fields[4]  # ALT allele
        qual = fields[5]  # QUAL
        filter = fields[6]  # FILTER
        info = fields[7]  # INFO
        lifted_info = info # will be updated if needed

        var_head = "\t".join(fields[:7])

        svtype = normalize_svtype(chrom, pos, ref, alt)

        # Optional FORMAT and sample columns
        ####################################
        # 'FORMAT' column from VCF is stored in 'fmt' to avoid conflict
        # with Python's built-in 'format' function  
        fmt = fields[8] if len(fields) > 8 else ""  # FORMAT field
        samples = fields[9:] if len(fields) > 9 else []  # Sample columns

        return cls(
            chrom=chrom,
            pos=pos,
            id=raw_id,
            lifted_id=id,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=filter,
            info=info,
            format=fmt,
            samples=samples,
            line_number=line_number,
            var_head=var_head,
            svtype=svtype,
            lifted_info=lifted_info
        )

    def to_vcf_line(self) -> str:
        """
        Convert the Variant object back into a VCF-formatted line.

        Returns:
            str: Tab-separated string representing the VCF line

        Notes:
            - Includes FORMAT and sample fields if present
            - Ensures proper order: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples...
        """
        fields = [
            self.lifted_chrom,
            str(self.lifted_pos),
            self.lifted_id,
            self.lifted_ref,
            self.lifted_alt,
            self.lifted_qual,
            self.lifted_filter,
            self.lifted_info
        ]

        # Add FORMAT column if exists
        if self.lifted_format:
            fields.append(self.lifted_format)

        # Add sample columns if exists
        if self.samples:
            fields.extend(self.lifted_samples)

        # Join all fields with tab characters
        return "\t".join(fields)
