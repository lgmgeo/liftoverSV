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

import re
from io_tools.file_utils import print_flush as print

# =========================
# Normalize SV type
# =========================

# SVtype in which category: DUP? DEL? INV? INS? TRA? None?
##########################################################
# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
# - Type2: angle-bracketed SV notation:   alt="<INS>", "<DEL>", ...
# - Type3: squared-bracketed SV notation: alt="G]17:1584563]" or alt="G]chr17:1584563]"
def normalize_svtype(chrom, pos, ref, alt):
    """Determine the structural variant type: DEL, DUP, INS, INV, TRA, or None."""
    svtype = "None"

    # angle-bracketed notation
    if re.search(r"del|loss|<CN0>|<CN1>", alt, re.IGNORECASE):
        svtype = "DEL"
    elif re.search(r"dup|gain|MCNV", alt, re.IGNORECASE):
        svtype = "DUP"
    elif re.search(r"<CN(\d+)>", alt, re.IGNORECASE):
        i = int(re.search(r"<CN(\d+)>", alt, re.IGNORECASE).group(1))
        if i > 1:
            svtype = "DUP"
    elif re.search(r"inv", alt, re.IGNORECASE):
        svtype = "INV"
    elif re.search(r"ins|MEI|alu|line|sva", alt, re.IGNORECASE):  # "DEL_ALU" is set to "DEL", OK!
        svtype = "INS"
    elif re.search(r"TRA|TRN", alt):
        svtype = "TRA"
    else:
        # Breakend notation [square brackets]?
        pattern = r"([acgtnACGTN]*)(\[|\])([^\:]+):(\d+)(\[|\])([acgtnACGTN]*)"
        m = re.match(pattern, alt)
        if m:
            # square-bracketed notation ]...]
            baseLeft, bracketLeft, inBracketChrom, inBracketStart, bracketRight, baseRight = m.groups()
            inBracketStart = int(inBracketStart)
            if chrom != inBracketChrom:
                # TRA (chrom_#CHROM != chrom_ALT)
                #################################
                # Example:
                # 17      198982  trn_no_mateid_a A       A]2:321681]
                # 2       321681  trn_no_mateid_b G       G]17:198982]
                svtype = "TRA"
            elif len(baseLeft) > 1 or len(baseRight) > 1:
                # INS ("first mapped base" is followed by the inserted sequence)
                ################################################################
                # Example:
                # 13      53040041        ins_by_gridss   T       TATATATATACACAC[13:53040042[  => The "T" at the left corresponds to position 53040041
                # 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA  => The "A" at the right corresponds to position 53040042
                #                                                                               => The inserted sequence is "ATATATATACACAC"
                svtype = "INS"
            elif ((bracketRight == "]" and baseRight and inBracketStart > pos) or 
                  (bracketLeft == "[" and baseLeft and inBracketStart < pos)):
                # DUP ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; REF is outside of the brackets)
                #########################################################################################################
                # Example:
                # 2       3000    breakend_dup_a  T       ]2:5000]T
                # 2       5000    breakend_dup_b  T       T[2:3000[
                svtype = "DUP"
            elif ((bracketRight == "]" and baseLeft) or (bracketRight == "[" and baseRight)):
                    # INV ("first mapped base" is contained in the bracket: "N]" or "[N")
                #####################################################################
                # Example 1:
                # 3       2999    breakend_inv_1_a        T       T]3:5000]
                # 3       5000    breakend_inv_1_b        T       [3:2999[T
                # Example 2:
                # 3       3000    breakend_inv_2_a        T       [3:5001[T
                # 3       5001    breakend_inv_2_b        T       T]3:3000]
                # Example 3:
                # 3       3000    breakend_inv_3_a        T       [3:5001[T
                # 3       5001    breakend_inv_3_b        T       [3:3000[T
                svtype = "INV"
            elif ((bracketRight == "[" and baseLeft and inBracketStart > pos) or
                  (bracketRight == "]" and baseRight and inBracketStart < pos)):
                # DEL ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; "first mapped base" is before the brackets
                ####################################################################################################################
                # Example:
                # 12      3000    breakend_del_1_a        T       T[12:5000[
                # 12      5000    breakend_del_1_b        T       ]12:3000]T
                svtype = "DEL"
        elif re.fullmatch(r"[acgtnACGTN.*]*", ref+alt):
            # e.g.: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
        
            # The GRIDSS author says that a . followed by bases refers to a single breakend where the reads cannot be uniquely mapped back.
            # e.g.: 2       39564894        gridss28_45b    T       .TTCTCTCATAACAAACCATGACATCCAGTCATTTAATACAATATGTCTGGGGTGGCTGGGCCCCTTTTTT 246.24  LOW_QUAL        
            refbis = re.sub(r"[.*]", "", ref)
            altbis = re.sub(r"[.*]", "", alt)
            variant_length = len(altbis) - len(refbis)
            if variant_length > 0:
                # insertion
                svtype = "INS"
                if ref not in alt:
                    svtype = "None"
            else:
                # deletion
                svtype = "DEL"
                if alt not in ref:
                    svtype = "None"

    return svtype
