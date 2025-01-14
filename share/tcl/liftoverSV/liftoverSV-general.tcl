############################################################################################################
# LiftOver_SVs_from_VCF 0.1.1_beta                                                                         #
#                                                                                                          #
# Copyright (C) 2024-current Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             #
#                                                                                                          #
# This program is free software; you can redistribute it and/or                                            #
# modify it under the terms of the GNU General Public License                                              #
# as published by the Free Software Foundation; either version 3                                           #
# of the License, or (at your option) any later version.                                                   #
#                                                                                                          #
# This program is distributed in the hope that it will be useful,                                          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


proc ContentFromFile {{File ""}} {
    if {[string equal $File ""]} {return ""}
    set f     [open $File r]
    set Text [read -nonewline $f]
    close $f
    return $Text
}

proc LinesFromFile {{File ""}} {
    return [split [ContentFromFile $File] "\n"]
}

proc WriteTextInFile {texte fichier} {
    set    fifi [open $fichier a]
    puts  $fifi $texte
    close $fifi
    return 1
}

proc isAnEmptyFile {VCFfile} {

    # Return 1 if the extension is not .vcf
    if {![regexp -nocase "\\.(vcf(.gz)?)$" $VCFfile]} {return 1}

    # After filtering lines beginning with "#", at least 1 line should be still present
    if {[regexp -nocase ".gz$" $VCFfile]} {
        set f [open "|gzip -cd $VCFfile"]
    } else {
        set f [open "$VCFfile"]
    }
    set test 0
    while {![eof $f]} {
        set L [gets $f]
        if {$L eq ""} {continue}
        if {[regexp "^#" $L]} {continue}
        incr test
        break
    }

    if {$test} {return 0} else {return 1}
}


# SVtype in which category: DUP? DEL? INV? INS? TRA? None?
##########################################################
# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG" 
# - Type2: angle-bracketed SV notation:   alt="<INS>", "<DEL>", ...
# - Type3: squared-bracketed SV notation: alt="G]17:1584563]" or alt="G]chr17:1584563]" 
proc normalizeSVtype {chrom pos ref alt} {

	set SVtype "None"

	# angle-bracketed notation <...>
    if {[regexp -nocase "del|loss|<CN0>|<CN1>" $alt]} {
        set SVtype "DEL"
    } elseif {[regexp -nocase "dup|gain|MCNV" $alt ]} {
        set SVtype "DUP"
    } elseif {[regexp -nocase "<CN(\[0-9\]+)>" $alt match i]} {
        if {$i>1} {set SVtype "DUP"}
    } elseif {[regexp -nocase "inv" $alt]} {
        set SVtype "INV"
    } elseif {[regexp -nocase "ins|MEI|alu|line|sva" $alt]} { ;# "DEL_ALU" is set to "DEL", OK!
        set SVtype "INS"
    } elseif {[regexp -nocase "TRA|TRN" $alt]} {
        set SVtype "TRA"
    } elseif {[regexp "(\[acgtnACGTN\]*)(\\\[|\\\])(\[^:\]+):(\[0-9\]+)(\\\[|\\\])(\[acgtnACGTN\]*)" $alt match baseLeft bracketLeft inBracketChrom inBracketStart bracketRight baseRight]} {
        # square-bracketed notation ]...]
        if {$chrom != $inBracketChrom} {
	        # TRA (chrom_#CHROM != chrom_ALT)
            #################################
            # Example:
            # 17      198982  trn_no_mateid_a A       A]2:321681]  
            # 2       321681  trn_no_mateid_b G       G]17:198982] 
			set SVtype "TRA"
		} elseif {[string length $baseLeft] > 1 || [string length $baseRight] > 1} {
            # INS ("first mapped base" is followed by the inserted sequence)
            ################################################################
            # Example:
            # 13      53040041        ins_by_gridss   T       TATATATATACACAC[13:53040042[  => The "T" at the left corresponds to position 53040041
            # 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA  => The "A" at the right corresponds to position 53040042
            #                                                                               => The inserted sequence is "ATATATATACACAC""
			set SVtype "INS"
        } elseif {($bracketRight eq "]" && $baseRight ne "" && $inBracketStart > $pos) || ($bracketLeft eq "\[" && $baseLeft ne "" && $inBracketStart < $pos)} {
            # DUP ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; REF is outside of the brackets)
            #########################################################################################################
            # Example:
            # 2       3000    breakend_dup_a  T       ]2:5000]T 
            # 2       5000    breakend_dup_b  T       T[2:3000[     
			set SVtype "DUP"
        } elseif {($bracketRight eq "]" && $baseLeft ne "") || ($bracketRight eq "\[" && $baseRight ne "")} {
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
			set SVtype "INV"
        } elseif {($bracketRight eq "\[" && $baseLeft ne "" && $inBracketStart > $pos) || ($bracketRight eq "]" && $baseRight ne "" && $inBracketStart < $pos)} {
            # DEL ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; "first mapped base" is before the brackets
            ####################################################################################################################
            # Example:
            # 12      3000    breakend_del_1_a        T       T[12:5000[   
            # 12      5000    breakend_del_1_b        T       ]12:3000]T 
			set SVtype "DEL"
		}
    } elseif {[regexp -nocase "^\[acgtnACGTN.*\]+$" $ref$alt]} {
		# e.g.: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"

        # The GRIDSS author says that a . followed by bases refers to a single breakend where the reads cannot be uniquely mapped back.
        # e.g.: 2       39564894        gridss28_45b    T       .TTCTCTCATAACAAACCATGACATCCAGTCATTTAATACAATATGTCTGGGGTGGCTGGGCCCCTTTTTT 246.24  LOW_QUAL

        regsub -all "\[*.\]" $ref "" refbis 
        regsub -all "\[*.\]" $alt "" altbis

        set variantLength [expr {[string length $altbis]-[string length $refbis]}]

        if {$variantLength>0} {
            # insertion
            set SVtype "INS"
        } else {
            # deletion
            set SVtype "DEL"
        }

        #if {[regexp "\\\." $ref] || [regexp "\\\." $alt]} {
		#	# GRIDSS BAD QUALITY          
		#} 
 
        if {[string length $refbis] < 10000 && [string length $altbis] < 10000} { ;# else we can have "couldn't compile regular expression pattern: out of memory" in the next regexp!!!
            if {![regexp -nocase "$refbis" $altbis] && ![regexp -nocase "$altbis" $refbis]} {
               # Complex SV: AGT>ATTGCATGGACCTGAGTCCCCAAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCGGGGGGGGGG
               set SVtype "None"
            }
        } 
	}


    return $SVtype
}

# Check for a VCF file (.vcf or .vcf.gz) or for a chain file
############################################################
# => Return "with" for a file with a "chr" prefix
# => Return "without" for a file without a "chr" prefix
proc fileWithChr {fileToCheck} {

	if {[regexp -nocase "vcf.gz$" $fileToCheck]} {
		set f [open "| gzip -cd $fileToCheck"]
		set type "vcf"
    } elseif {[regexp -nocase "vcf$" $fileToCheck]} {
        set f [open "$fileToCheck"]
        set type "vcf"
    } elseif {[regexp -nocase ".chain$" $fileToCheck]} {
        set f [open "$fileToCheck"]
        set type "chain"
    }

	if {$type eq "vcf"} {
	    while {![eof $f]} {
	        set L [gets $f]
			if {[regexp "^#" $L]} {continue}
			if {[string range [lindex $L 0]	0 2] eq "chr"} {
				return "with"
			} else {
				return "without"
			}
		}
	} elseif {$type eq "chain"} {
        while {![eof $f]} {
            set L [gets $f]
            if {[regexp "^#" $L]} {continue}
			if {[regexp "chr" $L]} {
                return "with"
            } else {
                return "without"
            }
		}
	}
}


