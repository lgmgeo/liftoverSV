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

import re
from typing import Dict, Any, Tuple, Optional
from io_tools.file_utils import print_flush as print
from io_tools.chain_lifter import ChainLifter
from io_tools.fasta_extractor import FastaExtractor
from core.variant import Variant  


class LiftoverEngine:
    """
    Engine to lift over VCF variants from a query genome to a target genome.

    Uses ChainLifter for coordinate conversion and FastaExtractor
    for sequence extraction. Handles drop cases and outputs
    mapped and unmapped variants.
    """

    def __init__(self, chain: ChainLifter, extractor: FastaExtractor, percent: int):
        self.chain = chain
        self.extractor = extractor
        self.percent = percent

        # Counters for unmapped statistics
        self.n_mapped = 0
        self.n_unmapped = 0
        self.case_counts = {f"case{i}": 0 for i in range(1, 6)}

        # Contig lines from the header: 
        # ##contig=<ID=chr1,length=...
        # => will be handle in "add_new_header_lines_and_sort"      

        # Memorize all the contigs in the VCF after the lift
        self.S_lifted_contigs = set()

        # Memorize if there are new INFO/FORMAT/FILTER annotations to add in the header lines
        self.S_SVlines_INFO, self.S_SVlines_FORMAT, self.S_SVlines_FILTER = set(), set(), set()


    def lift_variant(self, v: Variant, g_liftoverSV: Dict[str, Any]) -> Tuple[Optional[str], Optional[str]]:
        """
        Lift a single Variant object. 

        Returns: 1) a "lifted VCF line" or None
                 2) a reason for unmapped cases or None
        """

        # Lift POS
        new_chrom, new_pos = self.chain.lift(v.chrom, v.pos)
        if new_chrom is None:
            self.n_unmapped += 1
            self.case_counts["case1"] += 1
            return None, f"{v.var_head}\tPOS not lifted"
        v.lifted_chrom = new_chrom
        v.lifted_pos = new_pos
        
        #######################################
        # Lift over REF and ALT:
        # - REF contains '.' or '*' inside the sequence (not at the pos or end)?
        # - Coordinates of the last NT in REF
        # - check sequence ALT notation
        # - check square-bracketed notation
        # - check angle-bracketed notation
        #######################################

        # VCF format: The REF and ALT Strings must include the base before the variant (which must be reflected in the POS field), unless the variant occurs
        #             at position 1 on the contig in which case it must include the base after the variant

        # Badly formatted REF?
        if "." in v.ref[1:-1] or "*" in v.ref[1:-1]:
            # REF contains '.' or '*' inside the sequence (not at the pos or end)
            # => not lifted, skip
            # NOTE: if ref="." then ref[1:-1]=""
            self.n_unmapped += 1
            self.case_counts["case5"] += 1
            return None, f"{v.var_head}\tREF badly formatted (contains '.' or '*' inside the sequence)"
        
        # Lift REF: Coordinates of the last NT
        ref_clean = re.sub(r"[.*]", "", v.ref)
        v.ref_length = len(ref_clean)
        v.last_nt_coord = v.pos + v.ref_length - 1
        new_last_nt_chrom, new_last_nt_coord = self.chain.lift(v.chrom, v.last_nt_coord)
        if new_last_nt_coord is None:
            self.n_unmapped += 1
            self.case_counts["case1"] += 1
            return None, f"{v.var_head}\tLast NT of REF not lifted"
        if new_last_nt_chrom != new_chrom:
            self.n_unmapped += 1
            self.case_counts["case2"] += 1
            return None, f"{v.var_head}\tPOS and last NT of REF are lifted on different chromosomes ({new_chrom}, {new_last_nt_chrom})"
        v.lifted_last_nt_chrom = new_last_nt_chrom
        v.lifted_last_nt_coord = new_last_nt_coord

        if re.match(r"^[acgtnACGTN.*]+$", v.ref + v.alt):
            # Sequence ALT notation
            reason = self._lift_sequence_alt_notation(v)
            if reason:
                return None, reason
        else:
            # Square-bracketed or angle-bracketed ALT notation
            reason = self._lift_bracketed_alt_notation(v)
            if reason:  
                return None, reason
            
            # Lift over the square-bracketed ALT
            match = re.match(r"([ACGTNacgtn.*]*)(\[|\])([^:]+):([0-9]+)(\[|\])([AC.*GTNacgtn.*]*)", v.alt)
            if match:
                # Square-bracketed ALT notation
                reason = self._lift_square_bracketed_alt_notation(v, match)
                if reason:  
                    return None, reason



        # Lift over END/SVEND
        #####################
        end = -10
        match = re.search(r'(^END|;END|^SVEND|;SVEND)=(\d+)(;|$)', v.lifted_info)
        if match:
            end = int(match.group(2))
            new_chrom_end, new_end = self.chain.lift(v.chrom, end)
            if new_end is None:
                # Case1: END not lifted, skip
                self.n_unmapped += 1
                self.case_counts["case1"] += 1
                return None, f"{v.var_head}\tEND ({v.chrom}:{end}) not lifted"
                                
            if v.lifted_chrom != new_chrom_end and v.svtype not in ("TRA", "None"):
                # Case 2: lifted POS / ALT / END are on different chromosomes
                self.n_unmapped += 1
                self.case_counts["case2"] += 1
                return None, f"{v.var_head}\tPOS and END are lifted on different chrom ({v.lifted_chrom} # {new_chrom_end})"

            if new_end < v.lifted_pos:
                # Case3: reversed coordinates order between POS and END during the lift
                self.n_unmapped += 1
                self.case_counts["case3"] += 1
                return None, f"{v.var_head}\tlifted_POS ({v.lifted_pos}) > lifted_END ({new_end})"

            svlen = end - v.pos
            svlen_lifted = new_end - v.lifted_pos
            if svlen_lifted < svlen * (1 - self.percent) or svlen_lifted > svlen * (1 + self.percent):
                # Case4: The distance between the two lifted positions changes significantly (Default: difference between both SVLENs > 5%)
                self.n_unmapped += 1
                self.case_counts["case4"] += 1
                return None, f"{v.var_head}\tthe distance between lifted_END ({new_end}) and lifted_POS ({v.lifted_pos}) changes significantly (svlen diff > {self.percent} %)"               

            # Update info with the new END coordinate
            info = re.sub(r'(^END|;END)=([^;]+)(;|$)', f'\\1={new_end}\\3', v.lifted_info)
            v.lifted_info = info

        # Lift over INFO/SVLEN, INFO/SVSIZE
        ###################################
        # Set SVLEN/SVSIZE to "." for SV type not in "DUP", "DEL", "INV" or "INS" (=> TRA, CPX...)
        if v.svtype not in {"DUP", "DEL", "INV", "INS"}: 
            # not DEL, DUP, INV or INS
            if 'svlen_lifted' in locals():
                svlen_lifted = "."

        # Keep the same SVLEN/SVSIZE for INS (the number of the inserted bases remains the same)
        if v.svtype == "INS":
            if 'svlen_lifted' in locals():
                del svlen_lifted

        # Lift over INFO/SVLEN, INFO/SVSIZE
        if 'svlen_lifted' in locals():
            info = re.sub(r'(^SVLEN|;SVLEN)=(-)?(\d+)(;|$)', lambda m: f"{m.group(1)}={m.group(2) or ''}{svlen_lifted}{m.group(4)}", v.lifted_info)
            v.lifted_info = info
            info = re.sub(r'(^SVSIZE|;SVSIZE)=(-)?(\d+)(;|$)', lambda m: f"{m.group(1)}={m.group(2) or ''}{svlen_lifted}{m.group(4)}", v.lifted_info)
            v.lifted_info = info

        # Clean up
        if 'svlen_lifted' in locals():
            del svlen_lifted
        
        # Lift over INFO/CIPOS, INFO/CIEND:
        ###################################
        #
        # If present, the number of entries must be twice the number of ALT alleles. CIEND consists of successive pairs
        # of records encoding the confidence interval start and end offsets relative to the END position inferred by SVLEN
        # for each ALT allele. For symbolic structural variants, the first in the pair must not be greater than 0, and the second must not be less than 0.
        # => "$cipos_1 <= 0"  and "$cipos_2 >= 0"
        #
        # Check and modify if needed the CIPOS/CIEND values in order to have:
        #   => POS-CIPOS > 0
        #   => END+CIEND < chrom_length

        # CIPOS
        m = re.search(r'(^CIPOS|;CIPOS)=([^;]+)(;|$)', v.lifted_info)
        if m:
            cipos_vals = [int(x) for x in m.group(2).split(",")]
            cipos_1, cipos_2 = cipos_vals[0], cipos_vals[1]

            if cipos_1 <= 0 and cipos_2 >= 0:
                if v.lifted_pos + cipos_1 <= 0:
                    cipos_1 = -v.lifted_pos + 1
                    # Replace CIPOS in info
                    info = re.sub(
                        r'(^CIPOS|;CIPOS)=([^;]+)(;|$)',
                        f"{m.group(1)}={cipos_1},{cipos_2}{m.group(3)}",
                        v.lifted_info
                    )
                    v.lifted_info = info

        # CIEND
        m = re.search(r'(^CIEND|;CIEND)=([^;]+)(;|$)', v.lifted_info)
        if m:
            ciend_vals = [int(x) for x in m.group(2).split(",")]
            ciend_1, ciend_2 = ciend_vals[0], ciend_vals[1]

            if ciend_1 <= 0 and ciend_2 >= 0:
                chrom_size = g_liftoverSV['size_chrom_target'][v.lifted_chrom]
                if new_end + ciend_2 > chrom_size:
                    ciend_2 = chrom_size - new_end
                    # Replace CIEND in info
                    info = re.sub(
                        r'(^CIEND|;CIEND)=([^;]+)(;|$)',
                        f"{m.group(1)}={ciend_1},{ciend_2}{m.group(3)}",
                        v.lifted_info
                    )
                    v.lifted_info = info

        # Update S_lifted_contigs
        self.S_lifted_contigs.add(v.lifted_chrom)

        # Memorize all the FILTER annotations presents in the SV lines
        self.S_SVlines_FILTER.update(v.filter.split(";"))

        # Memorize all the INFO annotations presents in the SV lines
        infos_to_check = re.sub(r"=[^;]+", "", v.lifted_info)  
        self.S_SVlines_INFO.update(infos_to_check.split(";"))
    
        # Memorize all the FORMAT annotations presents in the SV lines
        if v.format:  # only if not empty, else add {''}
            self.S_SVlines_FORMAT.update(set(v.format.split(":")))     

        # Return a lifted VCF line
        self.n_mapped += 1
        fields = [
            v.lifted_chrom,
            v.lifted_pos,
            v.lifted_id,
            v.lifted_ref,
            v.lifted_alt,
            v.qual,
            v.filter,
            v.lifted_info,
            v.format,
            *v.samples
        ]
        vcf_line = "\t".join(map(str, fields))

        return vcf_line, None 
    



    def _lift_sequence_alt_notation(self, v: Variant) -> Optional[str]:
        """
        Handle sequence ALT notation (ACGTN, . et *)

        Update:
            v.lifted_ref
            v.lifted_alt
        Returns:
            None if the ALT notation is valid.
            A string reason if unmapped.
        """
        
        if "." in v.alt[1:-1] or "*" in v.alt[1:-1]:
            # ALT contains '.' or '*' inside the sequence (not at the pos or end)
            # => not lifted       
            self.n_unmapped += 1
            self.case_counts["case5"] += 1
            return f"{v.var_head}\tALT badly formatted (contains '.' or '*' inside the sequence)"
   
        alt_clean = re.sub(r"[.*]", "", v.alt)
        alt_length = len(alt_clean)

        if re.match(r"^N+$", v.ref) or v.ref == ".":
            # REF is composed only of "N" or REF = "."
            v.lifted_ref = v.ref
            v.lifted_alt = v.alt

        elif re.match(r"^N+$", v.alt) or v.alt == ".":
            # ALT is composed only of "N" or ALT = "."
            v.lifted_alt = v.alt
            new_reftmp = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
            # Keep "." or "*" if any (though not expected in REF)
            # (=> replace the sequence from "ref" with the lifted sequence)
            v.lifted_ref = re.sub(r"[^.*]+", new_reftmp, v.ref)

        else:
            # REF and ALT are different of "." and not composed only of "N"
            if v.ref_length < alt_length:
                # Insertion

                # "re.escape(ref)"" escapes all special regex characters in the string 'ref',
                # so that 'ref' is treated as a literal string in the regular expression.
                # => characters like '.', '*', '[', ']' lose their special meaning.
                # "re.match" always tries to match the pattern at the beginning of the string.
                m = re.match(re.escape(v.ref), v.alt, flags=re.IGNORECASE)
                if m:
                    # Insertion where REF is at the beginning of ALT
                    # e.g. REF="T" and ALT="TATTCTTG"
                    insertion = v.alt[m.end():]
                    new_reftmp = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
                    v.lifted_ref = re.sub(r"[^.*]+", new_reftmp, v.ref)
                    v.lifted_alt = f"{new_reftmp}{insertion}"

                elif v.alt.startswith(".") and v.alt.endswith(v.ref):
                    # Insertion with a Single Breakend: 
                    # the REF sequence is at the opposite side of the "." in the ALT sequence 
                    # e.g. REF="G" and ALT=".TTTTG"
                    insertion = v.alt[:-len(v.ref)]
                    v.lifted_ref = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
                    v.lifted_alt = f"{insertion}{v.lifted_ref}"

                else:
                    self.n_unmapped += 1
                    self.case_counts["case5"] += 1
                    return f"{v.var_head}\tcomplex sequence notation. INS: REF not at the beginning of ALT"
                    

            elif v.ref_length > alt_length:
                # Deletion
                m = re.match(re.escape(v.alt), v.ref, flags=re.IGNORECASE)
                if m:
                    # Deletion where ALT is at the beginning of REF
                    # e.g. REF="TATTCTTG" and ALT="T"
                    deletion = v.ref[m.end():]
                    last_nt_before_deletion = v.pos + v.ref_length - len(deletion) - 1
                    new_last_nt_before_deletion_chrom, new_last_nt_before_deletion = self.chain.lift(v.chrom, last_nt_before_deletion)
                    if new_last_nt_before_deletion is None:
                        self.n_unmapped += 1
                        self.case_counts["case1"] += 1
                        return f"{v.var_head}\tLast NT before deletion not lifted (see REF and ALT) (DEL)"

                    if new_last_nt_before_deletion_chrom != v.lifted_chrom:
                        self.n_unmapped += 1
                        self.case_counts["case2"] += 1
                        return f"{v.var_head}\tlifted_#CHROM ({new_chrom}) and lifted_last_NT_before_del ({new_last_nt_before_deletion_chrom}) are located on different chromosomes"

                    new_alt_tmp = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, new_last_nt_before_deletion)
                    v.lifted_alt = re.sub(r"[^.*]+", new_alt_tmp, v.alt)
                    v.lifted_ref = f"{v.lifted_alt}{deletion}"

                else:
                    self.n_unmapped += 1
                    self.case_counts["case5"] += 1
                    return f"{v.var_head}\tcomplex sequence notation. DEL: ALT not at the beginning of REF"

            else:
                # Equal lengths -> possible breakend or substitution
                if f".{v.ref}" == v.alt:
                    # Single breakend (e.g. REF = A and ALT = .A)
                    v.lifted_ref = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
                    v.lifted_alt = f".{v.lifted_ref}"

                elif f"{v.ref}." == v.alt:
                    # Single breakend (e.g. REF = G and ALT = G.)
                    v.lifted_ref = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
                    v.lifted_alt = f"{v.lifted_ref}."

                else:
                    # Complex sequence notation or substitution

                    # Check substitution (<= 50 bp): 
                    if len(v.ref) <= 50:
                        v.lifted_ref = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, v.lifted_last_nt_coord)
                        if v.ref != v.lifted_ref:
                            # The REF sequence differs from the original after liftover.
                            self.n_unmapped += 1
                            self.case_counts["case5"] += 1
                            return f"{v.var_head}\tthe REF sequence differs from the original after liftover"
                        else:
                            # The REF sequence is the same than the original after liftover
                            v.lifted_alt = v.alt
                
        # Variant lifted successfully
        return None
                    

    def _lift_bracketed_alt_notation(self, v: Variant) -> Optional[str]:
        """
        - Square-bracketed notation ([chr:pos[ ou ]chr:pos])
        - Angle-bracketed notation (<SVTYPE>)
        """
        # Do not lift over ALT (e.g. REF=A and ALT=<DEL>)
        # Initially, keep ALT as is
        v.lifted_alt = v.alt ;# it will be updated later for square-bracketed ALT notation

        # Lift over the REF sequence (e.g. REF=A and ALT=<DEL>)
        if re.fullmatch(r"N+", v.ref):
            # REF is composed entirely of 'N'
            v.lifted_ref = v.ref
        elif v.ref == ".":
            # REF = "."
            v.lifted_ref = "."
        else:
            # Extract the sequence from the lifted genome coordinates
            v.lifted_ref = self.extractor.get_sequence(v.lifted_chrom, v.lifted_pos, int(v.lifted_last_nt_coord))

    def _lift_square_bracketed_alt_notation(self, v: Variant, match: re.Match) -> Optional[str]:
        """
        Square-bracketed ALT notation ([chr:pos[ ou ]chr:pos])
        """
        # match = re.match(r"([ACGTNacgtn.*]*)(\[|\])([^:]+):([0-9]+)(\[|\])([ACGTNacgtn.*]*)", v.alt)
        base_left, bracket_left, chrom_alt, pos_alt, bracket_right, base_right = match.groups()
        pos_alt = int(pos_alt)
        
        new_chrom_alt, new_pos_alt = self.chain.lift(chrom_alt, pos_alt)
        if new_pos_alt is None:
            # Case1: not lifted
            self.n_unmapped += 1
            self.case_counts["case1"] += 1
            return f"{v.var_head}\tALT not lifted"

        if v.svtype not in ("TRA", "None"):
            # Case 2: lifted POS / ALT are on different chromosomes
            if new_chrom_alt != v.lifted_chrom:
                self.n_unmapped += 1
                self.case_counts["case2"] += 1
                return f"{v.var_head}\tifted_#CHROM ({v.lifted_chrom}) and lifted_alt_chrom ({new_chrom_alt}) are located on different chromosomes"
                

            # Case3: reversed coordinates order between POS and ALT
            if (v.lifted_pos < new_pos_alt and v.pos > pos_alt) or (v.lifted_pos > new_pos_alt and v.pos < pos_alt):
                self.n_unmapped += 1
                self.case_counts["case3"] += 1
                return f"{v.var_head}\tlifted_POS ({v.lifted_pos}) and lifted_ALT ({new_pos_alt}) are in reverse order relative to POS ({v.pos}) and ALT ({pos_alt})"
                

            # Case4: The distance between the two lifted positions changes significantly (Default: difference between both SVLENs > 5%)
            svlen = abs(pos_alt - v.pos)
            svlen_lifted = abs(new_pos_alt - v.lifted_pos)
            if svlen_lifted < svlen * (1 - self.percent) or svlen_lifted > svlen * (1 + self.percent):
                self.n_unmapped += 1
                self.case_counts["case4"] += 1
                return f"{v.var_head}\tthe distance between lifted_ALT ({new_pos_alt}) and lifted_POS ({v.lifted_pos}) changes significantly (svlen diff > {self.percent*100}%)"
                
    
        # Lift the sequence inside the square-bracketed ALT 
        # Examples:
        # A]chr2:32156] → "A" = v.lifted_ref
        # [chr2:32156[A → "A" = v.lifted_ref
        # Insertion example:
        # ACCCCC[chr2:32156[ → chr2:32156 corresponds to "A" (first base before the inserted sequence)
        # ]chr2:32156]CCCCCA → chr2:32156 corresponds to "A" (last base after the inserted sequence)
        if base_left and not base_right:
            # ALT of type base[chr:pos[
            if base_left == ".":
                new_base_left = "."
            else:
                new_base_left = v.lifted_ref
            v.lifted_alt = f"{new_base_left}{base_left[1:]}{bracket_left}{new_chrom_alt}:{new_pos_alt}{bracket_right}"

        elif base_right and not base_left:
            # ALT of type ]chr:pos]base
            if base_right == ".":
                new_base_right = "."
            else:
                new_base_right = v.lifted_ref
            v.lifted_alt = f"{bracket_left}{new_chrom_alt}:{new_pos_alt}{bracket_right}{base_right[:-1]}{new_base_right}"

        else:
            # Badly formatted ALT
            self.n_unmapped += 1
            self.case_counts["case5"] += 1
            return f"{v.var_head}\tsquare-bracketed ALT notation ({v.alt}) badly formatted."

        # Update S_lifted_contigs
        self.S_lifted_contigs.add(new_chrom_alt)
      
    

      