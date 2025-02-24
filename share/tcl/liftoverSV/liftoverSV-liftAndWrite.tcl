############################################################################################################
# liftoverSV 0.2.0_beta                                                                                    #
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


# Extract the different genomic coordinates from the input VCF (1-based) to create a BED file (0-based)
proc extractAllTheGenomicCoordinates {} {
    
    global g_liftoverSV
    global g_L_chrom
    
    
    puts "\n...reading [file tail $g_liftoverSV(INPUTFILE)] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    set tmpBEDfile "$g_liftoverSV(OUTPUTDIR)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"].tmp.bed"
    
    set L_toWrite {}
    set n_toWrite 0 ;#Do not use $i, this could cause a bug if all the coordinates of the SV line are already memorized (during the intermediate writing)
    set i 0
    set nSV 0
    set g_L_chrom(beforeLifted) {}
    
    
    
    if {[regexp ".gz$" $g_liftoverSV(INPUTFILE)]} {
        set F [open "|gzip -cd $g_liftoverSV(INPUTFILE)"]
    } else {
        set F [open "$g_liftoverSV(INPUTFILE)"]
    }
    while {[gets $F L]>=0} {
        
        if {[regexp "^#" $L]} {
            # Memorize the contigs defined in the header lines:
            # e.g. ##contig=<ID=chr1,length=249250621>
            if {[regexp "^##contig=<ID=(\[^,>\]+)" $L match contig]} {
                if {[lsearch -exac $g_L_chrom(beforeLifted) $contig] eq -1} {
                    lappend g_L_chrom(beforeLifted) $contig
                }
            }
            continue
        }
        
        set Ls [split $L "\t"]
        
        set chrom [lindex $Ls 0]
        set VCFstart [lindex $Ls 1]
        set BEDstart [expr {$VCFstart-1}] ;# 1-based => 0-based
        set ref [lindex $Ls 3]
        set alt [lindex $Ls 4]
        set infos [lindex $Ls 7]
        
        incr nSV
        
        # POS
        if {![info exists localMemory($chrom,$BEDstart)]} {
            set localMemory($chrom,$BEDstart) 1
            lappend L_toWrite "$chrom\t$BEDstart\t$VCFstart\t$i"
            incr i
            incr n_toWrite
        }
        
        # Coordinates of the last NT in the REF
        regsub -all "\[*.\]" $ref "" refbis
        set REFlength [string length $refbis]
        set lastNTfromREF [expr {$BEDstart+$REFlength-1}]
        if {![info exists localMemory($chrom,$lastNTfromREF)]} {
            set localMemory($chrom,$lastNTfromREF) 1
            lappend L_toWrite "$chrom\t$lastNTfromREF\t[expr {$lastNTfromREF+1}]\t$i"
            incr i
            incr n_toWrite
        }
        
        # For a deletion indicated in REF/ALT with a sequence notation (e.g. ATGCCC > ATG)
        # => coordinates of the last NT before the deletion (in our example, coordinates of the "G")
        if {[regexp "^\[acgtnACGTN.*\]+$" $ref$alt]} {
            regsub -all "\[*.\]" $alt "" altbis
            set REFlength [string length $refbis]
            set ALTlength [string length $altbis]
            
            if {![regexp "^N+$" $ref] && $ref ne "." && ![regexp "^N+$" $alt] && $alt ne "." } {
                # REF and ALT are different of "." and not composed only of "N"
                if {$REFlength > $ALTlength} {
                    # Deletion
                    if {[regsub -nocase "^$alt" $ref "" deletion]} {
                        set lastNTbeforeDeletion [expr {$BEDstart+$REFlength-[string length $deletion]-1}]
                        if {![info exists localMemory($chrom,$lastNTbeforeDeletion)]} {
                            set localMemory($chrom,$lastNTbeforeDeletion) 1
                            lappend L_toWrite "$chrom\t$lastNTbeforeDeletion\t[expr {$lastNTbeforeDeletion+1}]\t$i"
                            incr i
                            incr n_toWrite
                        }
                    }
                }
            }
        }
        
        # INFO/END
        if {[regexp "(^END|;END)=(\[^;\]+)(;|$)" $infos match tutu end titi]} {
            if {![info exists localMemory($chrom,[expr {$end-1}])]} {
                set localMemory($chrom,[expr {$end-1}]) 1
                lappend L_toWrite "$chrom\t[expr {$end-1}]\t$end\t$i"
                incr i
                incr n_toWrite
            }
        }
        
        # INFO/SVEND
        if {[regexp "(^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos match tutu svend titi]} {
            if {![info exists localMemory($chrom,[expr {$svend-1}])]} {
                set localMemory($chrom,[expr {$svend-1}]) 1
                lappend L_toWrite "$chrom\t[expr {$svend-1}]\t$svend\t$i"
                incr i
                incr n_toWrite
            }
        }
        
        # Square-bracketed ALT notation
        if {[regexp "(\\\[|\\\])(\[^:\]+):(\[0-9\]+)(\\\[|\\\])" $alt match tutu chromALT VCFposALT]} {
            set BEDposALT [expr {$VCFposALT-1}]
            if {![info exists localMemory($chromALT,$BEDposALT)]} {
                set localMemory($chromALT,$BEDposALT) 1
                lappend L_toWrite "$chromALT\t$BEDposALT\t$VCFposALT\t$i"
                incr i
                incr n_toWrite
            }
        }
        
        
        if {$n_toWrite eq 100000} {
            WriteTextInFile [join $L_toWrite "\n"] $tmpBEDfile
            set L_toWrite {}
            set n_toWrite 0
        }
    }
    close $F
    puts "\t* $nSV SV given in input"
    puts "\t* creation of [file tail $tmpBEDfile] with the different coordinates to lift (CHROM:POS, CHROM:END...)"
    puts "\t  ($i genomic coordinates to lift)"
    if {$L_toWrite ne ""} {
        WriteTextInFile [join $L_toWrite "\n"] $tmpBEDfile
    }
    
    return $tmpBEDfile
}



# Lift over these different genomic coordinates
proc liftoverGenomicCoordinates {tmpBEDfile} {
    
    global g_liftoverSV
    
    puts "\n...lifting over [file tail $tmpBEDfile] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    regsub ".bed" $tmpBEDfile ".lifted.bed" tmpBEDfileLifted
    set command "$g_liftoverSV(LIFTOVER) $tmpBEDfile $g_liftoverSV(CHAIN) $tmpBEDfileLifted $g_liftoverSV(OUTPUTDIR)/unMapped"
    
    if {[catch {eval exec $command} Message]} {
        puts "\t* command line:"
        puts "\t\t$command"
        puts "\t* output:"
        foreach L [split $Message "\n"] {
            puts "\t\t$L"
        }
    }
    
    
    file delete -force $g_liftoverSV(OUTPUTDIR)/unMapped
    return $tmpBEDfileLifted
}



# Memorize the lifted coordinate:s
# Read from a BED file (0-based) and memorize in a VCF format (1-based)
proc memorizeGenomicCoordinates {tmpBEDfile tmpBEDfileLifted} {
    
    global g_liftoverSV
    global g_ID
    global g_liftedCoord
    global g_L_chrom
    
    puts "\n...memorizing the lifted coordinates ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    set F [open "$tmpBEDfile"]
    while {[gets $F L]>=0} {
        set Ls [split $L "\t"]
        set chrom [lindex $Ls 0]
        
        set BEDpos [lindex $Ls 1]
        set VCFpos [expr {$BEDpos+1}]
        set id [lindex $Ls end]
        set g_ID($chrom,$VCFpos) $id
    }
    
    set g_L_chrom(lifted) {}
    set F [open "$tmpBEDfileLifted"]
    while {[gets $F L]>=0} {
        set Ls [split $L "\t"]
        set id [lindex $Ls end]
        set chrom [lindex $Ls 0]
        set BEDpos [lindex $Ls 1]
        set VCFpos [expr {$BEDpos+1}]
        set g_liftedCoord($id) "$chrom,$VCFpos"
        if {[lsearch -exac $g_L_chrom(lifted) $chrom] eq -1} {
            lappend g_L_chrom(lifted) $chrom
        }
    }
    file delete -force $tmpBEDfile
    file delete -force $tmpBEDfileLifted
}



# Write the lifted VCF
######################
# Rules:
# => lift #CHROM, POS, REF, ALT, INFO/END and INFO/SVEND
# => drop the SV if:
#   - Case1: one position (start or end) is lifted while the other doesn't
#   - Case2: one position (start or end) goes to a different chrom from the other (except for translocation)
#   - Case3: "lifted start" > "lifted end"
#   - Case4: The distance between the two lifted positions changes significantly (Default: difference between both SVLENs > 5%)
#   - Case5: The REF and ALT features are represented with complex sequences:
#		- Deletion: the REF sequence is not at the beginning of the ALT sequence e.g. `REF="ATTCTTG" and ALT="TC"
#		- Insertion: the ALT sequence is not at the beginning of the REF sequence e.g. `REF="TC" and ALT="ATTCTTG"
#       - Insertion with a Single Breakend: the REF sequence is not at the opposite side of the "." in the ALT sequence e.g. `REF="G" and ALT=".TTTTTTC"

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
#
# INFO/CIPOS and INF²O/CIEND:
############################
# Check and modify if needed the CIPOS/CIEND values in order to have:
#	=> POS-CIPOS > 0
#	=> END+CIEND < chrom_length
#
proc writeTheLiftedVCF {} {
    
    global g_liftoverSV
    global g_ID
    global g_liftedCoord
    global g_L_chrom
    
    
    if {![regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".unmapped" unmappedFile]} {
        puts "Check the output file name ($g_liftoverSV(OUTPUTFILE))."
        puts "Extension sould be \".vcf\""
        puts "Exit with error"
        exit 2
    }
    
    regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".tmp.vcf" tmpOUTPUTFILE
    puts "\n...writing temporary [file tail $tmpOUTPUTFILE] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    file delete -force $tmpOUTPUTFILE
    file delete -force $unmappedFile
    
    # SVLEN calculation only for the following svtypes (not INS, TRA, CPX...)
    set L_svtypes {DUP DEL INV}
    
    # Memorize if there are new INFO/FORMAT/FILTER annotations to add in the header lines
    set L_header_INFO ""
    set L_header_FORMAT ""
    set L_header_FILTER ""
    set L_SVlines_INFO ""
    set L_SVlines_FORMAT ""
    set L_SVlines_FILTER ""
    set L_new_INFO ""
    set L_new_FORMAT ""
    set L_new_FILTER ""
    
    set L_toWrite {}
    set at_least_1_SV_lifted 0
    set n_mapped 0
    set L_unmappedToWrite {}
    set n_unmapped 0
    
    set case1 0
    set case2 0
    set case3 0
    set case4 0
    set case5 0
    set case6 0
    
    if {[regexp ".gz$" $g_liftoverSV(INPUTFILE)]} {
        set F [open "| gzip -cd $g_liftoverSV(INPUTFILE)"]
    } else {
        set F [open "$g_liftoverSV(INPUTFILE)"]
    }
    
    set testContigs "toDo"
    
    while {[gets $F L]>=0} {
        
        # VCF header lines
        ##################
        if {[regexp "^#" $L]} {
            if {[regexp "^##contig=<ID=" $L]} {
                ##contig=<ID=chr1,length=249250621>
                regsub ",length=\[0-9\]+" $L "" L
                set testContigs "inProgress"
                lappend L_toWrite "$L"
                continue
            }
            if {$testContigs eq "inProgress"} {
                foreach chrom [lsort $g_L_chrom(lifted)] {
                    if {[lsearch -exact $g_L_chrom(beforeLifted) $chrom] eq -1} {
                        lappend L_toWrite "##contig=<ID=$chrom>"
                    }
                }
                set testContigs "done"
            }
            lappend L_toWrite "$L"
            
            ##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
            ##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
            ##FILTER=<ID=ID,Description="description">
            if {[regexp "^##INFO=<ID=(\[^,\]+)," $L match knownValue]} {
                lappend L_header_INFO $knownValue
            }
            if {[regexp "^##FORMAT=<ID=(\[^,\]+)," $L match knownValue]} {
                lappend L_header_FORMAT $knownValue
            }
            if {[regexp "^##FILTER=<ID=(\[^,\]+)," $L match knownValue]} {
                lappend L_header_FILTER $knownValue
            }
            
            # Rescue if there is no contig lines in the header
            if {[regexp "^#CHROM" $L] && $testContigs eq "toDo"} {
                foreach chrom [lsort $g_L_chrom(lifted)] {
                    if {[lsearch -exact $g_L_chrom(beforeLifted) $chrom] eq -1} {
                        lappend L_toWrite "##contig=<ID=$chrom>"
                    }
                }
                set testContigs "done"
            }
            
            continue
        }
        
        
        # SV lines
        ##########
        set Ls [split $L "\t"]
        
        set chrom [lindex $Ls 0]
        set start [lindex $Ls 1]
        set ID [lindex $Ls 2]
        set ref [lindex $Ls 3]
        set alt [lindex $Ls 4]
        set qual [lindex $Ls 5]
        set filter [lindex $Ls 6]
        set infos [lindex $Ls 7]
        set format_annot [lindex $Ls 8]
        
        set svtype [normalizeSVtype $chrom $start $ref $alt]
        # Drop SV if the SVTYPE is not DUP, DEL, INV, INS or TRA
        # if {$svtype eq "None"} {continue}
        # ...A voir si nécessaire
        
        
        # Lift over POS
        set theID $g_ID($chrom,$start)
        if {[info exists g_liftedCoord($theID)]} {
            set theNewChrom [lindex [split $g_liftedCoord($theID) ","] 0]
            set theNewStart [lindex [split $g_liftedCoord($theID) ","] 1]
        } else {
            # Case1
            lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tStart ($chrom,$start) not lifted"
            incr n_unmapped
            incr case1
            continue
        }
        
        #######################################
        # Lift over REF and ALT:
        # - square-bracketed notation
        # - angle-bracketed notation
        # - sequence notation
        #######################################
        
        # Case5: The REF and ALT features are represented with complex sequences:
        #   - Deletion: the REF sequence is not at the beginning of the ALT sequence e.g. `REF="ATTCTTG" and ALT="TC"
        #   - Insertion: the ALT sequence is not at the beginning of the REF sequence e.g. `REF="TC" and ALT="ATTCTTG"
        #   - Insertion with a Single Breakend: the REF sequence is not at the opposite side of the "." in the ALT sequence e.g. `REF="G" and ALT=".TTTTTTC"
        # => Drop SV with a complex sequence notation in REF/ALT
        
        # REF and ALT are represented with a sequence notation
        #
        # VCF format: The REF and ALT Strings must include the base before the variant (which must be reflected in the POS field), unless the variant occurs
        #             at position 1 on the contig in which case it must include the base after the variant
        if {[regexp "^\[acgtnACGTN.*\]+$" $ref$alt]} {
            
            # => sequence notation
            ######################
            
            regsub -all "\[*.\]" $ref "" refbis
            regsub -all "\[*.\]" $alt "" altbis
            set REFlength [string length $refbis]
            set ALTlength [string length $altbis]
            
            # Coordinates of the last NT in the REF
            set lastNTfromREF [expr {$start+$REFlength-1}]
            set theID $g_ID($chrom,$lastNTfromREF)
            set theNewlastNTfromREF [lindex [split $g_liftedCoord($theID) ","] 1]
            
            
            if {[regexp "^N+$" $ref] || $ref eq "."} {
                # REF is composed only of "N" or REF = "."
                set theNewREF $ref
                set theNewALT $alt
            } elseif {[regexp "^N+$" $alt] || $alt eq "."} {
                # ALT is composed only of "N" or ALT = "."
                set theNewALT $alt
                set theNewREFtmp [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
                regsub "\[^.*\]+" $ref "$theNewREFtmp" theNewREF ;# We keep the "." or "*" information (just in case. This should not be observed in REF, should only be observed in ALT)
            } else {
                # REF and ALT are different of "." and not composed only of "N"
                if {$REFlength < $ALTlength} {
                    # Insertion
                    if {[regsub -nocase "^$ref" $alt "" insertion]} {
                        set theNewREFtmp [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
                        regsub "\[^.*\]+" $ref "$theNewREFtmp" theNewREF
                        set theNewALT "$theNewREF$insertion"
                    } elseif {[regexp "^\\." $alt] && [regsub -nocase "$ref$" $alt "" insertion]} {
                        # Insertion with a Single Breakend: the REF sequence is at the opposite side of the "." in the ALT sequence e.g. REF="G" and ALT=".TTTTG"
                        set theNewREF [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
                        set theNewALT "$insertion$theNewREF"
                    } else {
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tcomplex sequence notation (see REF and ALT)"
                        incr n_unmapped
                        incr case5
                        continue
                    }
                } elseif {$REFlength > $ALTlength} {
                    # Deletion
                    if {[regsub -nocase "^$alt" $ref "" deletion]} {
                        set lastNTbeforeDeletion [expr {$start+$REFlength-[string length $deletion]-1}]
                        set theID $g_ID($chrom,$lastNTbeforeDeletion)
                        set theNewlastNTbeforeDeletion [lindex [split $g_liftedCoord($theID) ","] 1]
                        set theNewALTtmp [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTbeforeDeletion]
                        regsub "\[^.*\]+" $alt "$theNewALTtmp" theNewALT
                        set theNewREF "$theNewALT$deletion"
                    } else {
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tcomplex sequence notation (see REF and ALT)"
                        incr n_unmapped
                        incr case5
                        continue
                    }
                } else {
                    if {".$ref" eq $alt} {
                        # Single breakend (e.g. REF = A and ALT = .A)
                        set theNewREF [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
                        set theNewALT ".$theNewREF"
                    } elseif {"$ref." eq $alt} {
                        # Single breakend (e.g. REF = G and ALT = G.)
                        set theNewREF [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
                        set theNewALT "$theNewREF."
                    } else {
                        # Complex sequence notation or substitution
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tcomplex sequence notation (see REF and ALT)"
                        incr n_unmapped
                        incr case5
                        continue
                    }
                }
            }
        } else {
            
            # => square-bracketed and/or angle-bracketed ALT notation
            #########################################################
            
            # Do not lift over ALT (e.g. REF=A and ALT=<DEL>)
            set theNewALT $alt   ;# Will be erased after for square-bracketed ALT notation
            
            # Lift over the sequence in REF (e.g. REF=A and ALT=<DEL>)
            if {[regexp "^N+$" $ref]} {
                set theNewREF $ref       ;# REF is composed only of "N"
            } elseif {$ref eq "."} {
                set theNewREF "."        ;# REF = "."
            } else {
                # Coordinates of the last NT in the REF
                regsub -all "\[*.\]" $ref "" refbis
                set REFlength [string length $refbis]
                set lastNTfromREF [expr {$start+$REFlength-1}]
                set theID $g_ID($chrom,$lastNTfromREF)
                set theNewlastNTfromREF [lindex [split $g_liftedCoord($theID) ","] 1]
                set theNewREF [ExtractDNAseq $theNewChrom $theNewStart $theNewlastNTfromREF]
            }
        }
        
        
        # Memorize all the INFO annotations presents in the SV lines
        regsub -all "=\[^;\]+" $infos "" infosToCheck
        lappend L_SVlines_INFO {*}[split $infosToCheck ";"]
        set L_SVlines_INFO [lsort -unique $L_SVlines_INFO]
        
        # Memorize all the FORMAT annotations presents in the SV lines
        lappend L_SVlines_FORMAT {*}[split $format_annot ":"]
        set L_SVlines_FORMAT [lsort -unique $L_SVlines_FORMAT]
        
        # Memorize all the FILTER annotations presents in the SV lines
        lappend L_SVlines_FILTER {*}[split $filter ";"]
        set L_SVlines_FILTER [lsort -unique $L_SVlines_FILTER]
        
        # Lift over the square-bracketed ALT
        if {[regexp "(\[acgtnACGTN\]*)(\\\[|\\\])(\[^:\]+):(\[0-9\]+)(\\\[|\\\])(\[acgtnACGTN\]*)" $alt match baseLeft bracketLeft chromALT posALT bracketRight baseRight]} {
            set theIDalt $g_ID($chromALT,$posALT)
            if {[info exists g_liftedCoord($theIDalt)]} {
                set theNewChromALT [lindex [split $g_liftedCoord($theIDalt) ","] 0]
                set theNewPosALT [lindex [split $g_liftedCoord($theIDalt) ","] 1]
                if {$svtype ne "TRA" && $svtype ne "None"} {
                    if {$theNewChromALT ne $theNewChrom} {
                        # Case2
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_#CHROM ($theNewChrom) and lifted_ALT ($theNewChromALT) are located on different chromosomes"
                        incr n_unmapped
                        incr case2
                        continue
                    }
                    if {$theNewStart < $theNewPosALT && $start > $posALT} {
                        # Case3
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_POS ($theNewStart) < lifted_ALT ($theNewPosALT) while POS ($start) > ALT ($posALT)"
                        incr n_unmapped
                        incr case3
                        continue
                    }
                    if {$theNewStart > $theNewPosALT && $start < $posALT} {
                        # Case3
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_POS ($theNewStart) > lifted_ALT ($theNewPosALT) while POS ($start) < ALT ($posALT)"
                        incr n_unmapped
                        incr case3
                        continue
                    }
                    
                    set svlen [expr {abs($posALT-$start)}]
                    set svlenlifted [expr {abs($theNewPosALT-$theNewStart)}]
                    if {$svlenlifted < [expr {$svlen*(1-$g_liftoverSV(PERCENT))}] || $svlenlifted > [expr {$svlen*(1+$g_liftoverSV(PERCENT))}]} {
                        # Case4
                        lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe distance between lifted_ALT ($theNewPosALT) and lifted_POS ($theNewStart) changes significantly (svlen diff > $g_liftoverSV(PERCENT)%)"
                        incr n_unmapped
                        incr case4
                        continue
                    }
                }
            } else {
                # Case1
                lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tALT ($chromALT,$posALT) not lifted"
                incr n_unmapped
                incr case1
                continue
            }
            
            # Lift the sequence in the square-bracketed ALT
            # Examples:
            # A]chr2:32156] => "A" = $theNewREF
            # [chr2:32156[A => "A" = $theNewREF
            # INS:
            # ACCCCC[chr2:32156[ => chr2:32156 is "A" (first NT before the inserted sequence)
            # ]chr2:32156]CCCCCA => chr2:32156 is "A"(last NT after the inserted sequence)
            if {$baseLeft ne "" && $baseRight eq ""} {
                set theNewALT "$theNewREF[string range $baseLeft 1 end]$bracketLeft$theNewChromALT:$theNewPosALT$bracketRight"
            } elseif {$baseRight ne "" && $baseLeft eq ""} {
                set theNewALT "$bracketLeft$theNewChromALT:$theNewPosALT$bracketRight[string range $baseRight 0 end-1]$theNewREF"
            } else {
                lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe ALT features ($alt) is badly formatted."
                incr n_unmapped
                incr case5
                continue
            }
        }
        
        
        # Memorize the 6 first columns lifted
        set theNewL "$theNewChrom\t$theNewStart\t$ID\t$theNewREF\t$theNewALT\t$qual\t$filter"
        
        
        # Lift over END/SVEND
        #####################
        set end "-10"
        if {[regexp "(^END|;END|^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos match tutu end titi]} {
            set theID $g_ID($chrom,$end)
            if {[info exists g_liftedCoord($theID)]} {
                set theNewChromEND [lindex [split $g_liftedCoord($theID) ","] 0]
                set theNewEND [lindex [split $g_liftedCoord($theID) ","] 1]
                if {$theNewChrom ne $theNewChromEND} {
                    # Case2
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_start and lifted_END are located on different chrom ($theNewChrom # $theNewChromEND)"
                    incr n_unmapped
                    incr case2
                    continue
                }
                if {$theNewEND < $theNewStart} {
                    # Case3
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_start ($theNewStart) > lifted_end ($theNewEND)"
                    incr n_unmapped
                    incr case3
                    continue
                }
                
                set svlen [expr {$end-$start}]
                set svlenlifted [expr {$theNewEND-$theNewStart}]
                if {$svlenlifted < [expr {$svlen*(1-$g_liftoverSV(PERCENT))}] || $svlenlifted > [expr {$svlen*(1+$g_liftoverSV(PERCENT))}]} {
                    # Case4
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe distance between lifted_END ($theNewEND-) and lifted_POS ($theNewStart) changes significantly (svlen diff > $g_liftoverSV(PERCENT)%)"
                    incr n_unmapped
                    incr case4
                    continue
                }
            } else {
                # Case1
                lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tEND ($chrom,$end) not lifted"
                incr n_unmapped
                incr case1
                continue
            }
            
            regsub "(^END|;END)=(\[^;\]+)(;|$)" $infos "\\1=$theNewEND\\3" infos
        }
        
        
        # Lift over INFO/SVLEN, INFO/SVSIZE
        ###################################
        # Set SVLEN/SVSIZE to "." for SV type not equal to DEL, DUP, INV or INS (TRA, CPX...)
        if {[lsearch -exact $L_svtypes $svtype] eq -1} {
            if {$svtype ne "INS" && [info exists svlenlifted]} {set svlenlifted "."}
        }
        # Keep the same SVLEN/SVSIZE for insertion (the number of the inserted bases remains the same)
        if {$svtype eq "INS"} {catch {unset svlenlifted}}
        
        # Lift over INFO/SVLEN, INFO/SVSIZE for deletion, duplication and inversion
        if {[info exists svlenlifted]} {
            regsub "(^SVLEN|;SVLEN)=(-)?(\[0-9\]+)(;|$)" $infos "\\1=\\2$svlenlifted\\4" infos
            regsub "(^SVSIZE|;SVSIZE)=(-)?(\[0-9\]+)(;|$)" $infos "\\1=\\2$svlenlifted\\4" infos
        }
        catch {unset svlenlifted}
        
        
        # INFO/CIPOS and INFO/CIEND:
        ############################
        #
        # If present, the number of entries must be twice the number of ALT alleles. CIEND consists of successive pairs
        # of records encoding the confidence interval start and end offsets relative to the END position inferred by SV LEN
        # for each ALT allele. For symbolic structural variants, the first in the pair must not be greater than 0, and the second must not be less than 0.
        # => "$cipos_1 <= 0"  and "$cipos_2 >= 0"
        #
        # Check and modify if needed the CIPOS/CIEND values in order to have:
        #   => POS-CIPOS > 0
        #   => END+CIEND < chrom_length
        if {[regexp "(^CIPOS|;CIPOS)=(\[^;\]+)(;|$)" $infos match titi cipos]} {
            set cipos_1 [lindex [split $cipos ","] 0]
            # Lift over CIPOS only for good format
            if {$cipos_1 <= 0} {
                set cipos_2 [lindex [split $cipos ","] 1]
                if {$cipos_2 >= 0} {
                    if {[expr {$theNewStart+$cipos_1}] <= 0} {
                        set cipos_1 [expr {-$theNewStart+1}]
                        regsub "(^CIPOS|;CIPOS)=(\[^;\]+)(;|$)" $infos "\\1=$cipos_1,$cipos_2\\3" infos
                    }
                }
            }
        }
        if {[regexp "(^CIEND|;CIEND)=(\[^;\]+)(;|$)" $infos match titi ciend]} {
            set ciend_1 [lindex [split $ciend ","] 0]
            # Lift over CIEND only for good format
            if {$ciend_1 <= 0} {
                set ciend_2 [lindex [split $ciend ","] 1]
                if {$ciend_2 >= 0} {
                    if {[expr {$theNewEND+$ciend_2}] > $g_liftoverSV(sizeAfterLift,$theNewChrom)} {
                        set ciend_2 [expr {$g_liftoverSV(sizeAfterLift,$theNewChrom)-$theNewEND}]
                        regsub "(^CIEND|;CIEND)=(\[^;\]+)(;|$)" $infos "\\1=$ciend_1,$ciend_2\\3" infos
                    }
                }
            }
        }
        
        
        # tmp output (before the sort)
        ##############################
        append theNewL "\t$infos\t[join [lrange $Ls 8 end] "\t"]"
        
        lappend L_toWrite $theNewL
        set at_least_1_SV_lifted 1
        
        incr n_mapped
        
        if {![expr {$n_mapped%10000}]} {
            WriteTextInFile [join $L_toWrite "\n"] $tmpOUTPUTFILE
            set L_toWrite {}
        }
        if {$n_unmapped ne 0 && ![expr {$n_unmapped%10000}]} {
            WriteTextInFile [join $L_unmappedToWrite "\n"] $unmappedFile
            set L_unmappedToWrite {}
        }
        
    }
    close $F
    
    # Some SV are lifted
    if { $at_least_1_SV_lifted } {
        
        # Writings
        if {$L_toWrite ne ""} {
            WriteTextInFile [join $L_toWrite "\n"] $tmpOUTPUTFILE
        }
    }
    
    if {$L_unmappedToWrite ne ""} {
        WriteTextInFile [join $L_unmappedToWrite "\n"] $unmappedFile
    }
    puts "\t* $n_mapped mapped SV"
    
    if {$n_unmapped ne 0} {
        puts "\n...unmapped SV (see [file tail $unmappedFile] for details):"
        puts "\t* $n_unmapped unmapped SV"
        if {$case1} {puts "\t- $case1 SV with one position (start or end) lifted while the other doesn't"}
        if {$case2} {puts "\t- $case2 SV with one position (start or end) mapped to a different chrom from the other"}
        if {$case3} {puts "\t- $case3 SV with \"lifted start\" > \"lifted end\""}
        if {$case4} {puts "\t- $case4 SV with the distance between the two lifted positions that changes significantly (difference between both SVLENs > 5%)"}
        if {$case5} {puts "\t- $case5 SV with an ALT feature that is not a square/angle bracketed notation."}
    }
    
    # If no SV lifted (all are unmapped), stop here
    if {! $at_least_1_SV_lifted } {
        puts "\n...No SV lifted!"
        puts "\tExit without error"
        exit 0
    }
    
    # Add new header lines (only for INFO, FORMAT or FILTER) if needed
    foreach val $L_SVlines_INFO {
        if {[lsearch -exact $L_header_INFO $val] eq -1} {lappend L_new_INFO $val}
    }
    foreach val $L_SVlines_FORMAT {
        if {[lsearch -exact $L_header_FORMAT $val] eq -1} {lappend L_new_FORMAT $val}
    }
    foreach val $L_SVlines_FILTER {
        if {[lsearch -exact $L_header_FILTER $val] eq -1} {lappend L_new_FILTER $val}
    }
    
    
    # assembly / contig / FILTER / INFO / FORMAT
    addNewHeaderLinesAndSort $L_new_INFO $L_new_FORMAT $L_new_FILTER
    
    return
}


# Meta-information lines in VCF:
################################
##fileformat=VCFv4.4
##assembly=url
##contig=<ID=ctg1,length=81195210,URL=ftp://somewhere.example/assembly.fa,md5=f126cdf8a6e0c7f379d618ff66beb2da,...>
##FILTER=<ID=ID,Description="description">
##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
##ALT=<ID=type,Description="description">
##META=<ID=Tissue,Type=String,Number=.,Values=[Blood, Breast, Colon, Lung, ?]>
##SAMPLE=<ID=Sample1,Assay=WholeGenome,Ethnicity=AFR,Disease=None,Description="Patient germline genome from unaffected",DOI=url>
##PEDIGREE=<ID=TumourSample,Original=GermlineID>

# Add new header lines (only for INFO, FORMAT and FILTER) if needed:
####################################################################
# 1 - Retrieve information from $LIFTOVERSV/share/doc/liftoverSV/vcf_header_lines.txt
# 2 - If not present, as the format (Number, String) is not known, "Number=." and "Type=String" values are used by default:
#     - ##FORMAT=<ID=XXX,Number=.,Type=String,Description="XXX">
#     - ##INFO=<ID=YYY,Number=.,Type=String,Description="YYY">
#     - ##FILTER=<ID=YYY,Description="YYY">
proc addNewHeaderLinesAndSort {L_new_INFO L_new_FORMAT L_new_FILTER} {
    
    global g_liftoverSV
    
    
    puts "\n...checking the INFO, FORMAT and FILTER header lines ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    puts "\t* creation of $g_liftoverSV(OUTPUTFILE)"
    
    regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".tmp.vcf" tmpOUTPUTFILE
    file delete -force $g_liftoverSV(OUTPUTFILE)
    
    
    # Memorize information from $LIFTOVERSV/share/doc/liftoverSV/vcf_header_lines.txt
    set headerFile "$g_liftoverSV(installDir)/share/doc/liftoverSV/vcf_header_lines.txt"
    # Format in the headerFile:
    # - FILTER line: Key Description
    # - INFO line: Key Number Type Description
    # - FORMAT line: Field Number Type Description
    foreach L [LinesFromFile $headerFile] {
        if {[regexp "^#" $L]} {
            if {[regexp "^# FILTER"	$L]} {set feature "FILTER"}
            if {[regexp "^# INFO" $L]} {set feature "INFO"}
            if {[regexp "^# FORMAT" $L]} {set feature "FORMAT"}
            continue
        }
        if {$L eq ""} {continue}
        if {$feature eq "FILTER"} {
            set Key [lindex $L 0]
            set Description [ join [lrange $L 1 end] " "]
            set FILTER($Key,Description) $Description
        } elseif {$feature eq "INFO"} {
            set Key [lindex $L 0]
            set Number [lindex $L 1]
            set Type [lindex $L 2]
            set Description [ join [lrange $L 3 end] " "]
            set INFO($Key,Number) $Number
            set INFO($Key,Type) $Type
            set INFO($Key,Description) $Description
        } elseif {$feature eq "FORMAT"} {
            set Key [lindex $L 0]
            set Number [lindex $L 1]
            set Type [lindex $L 2]
            set Description [ join [lrange $L 3 end] " "]
            set FORMAT($Key,Number) $Number
            set FORMAT($Key,Type) $Type
            set FORMAT($Key,Description) $Description
        }
    }
    
    # Add INFO, FORMAT and FILTER
    foreach Key $L_new_INFO {
        if {[info exist INFO($Key,Description)]} {
            lappend L_INFO "##INFO=<ID=$Key,Number=$INFO($Key,Number),Type=$INFO($Key,Type),Description=\"$INFO($Key,Description)\">"
        } else {
            lappend L_INFO "##INFO=<ID=$Key,Number=.,Type=String,Description=\"$Key\">"
        }
    }
    foreach Key $L_new_FORMAT {
        if {[info exist FORMAT($Key,Description)]} {
            lappend L_FORMAT "##FORMAT=<ID=$Key,Number=$FORMAT($Key,Number),Type=$FORMAT($Key,Type),Description=\"$FORMAT($Key,Description)\">"
        } else {
            lappend L_FORMAT "##FORMAT=<ID=$Key,Number=.,Type=String,Description=\"$Key\">"
        }
    }
    foreach Key $L_new_FILTER {
        if {[info exist FILTER($Key,Description)]} {
            lappend L_FILTER "##FILTER=<ID=$Key,Description=\"$FILTER($Key,Description)\">"
        } else {
            lappend L_FILTER "##FILTER=<ID=$Key,Description=\"$Key\">"
        }
    }
    
    set L_fileformat ""
    set L_assembly ""
    set L_contig ""
    set L_ALT ""
    set L_META ""
    set L_SAMPLE ""
    set L_PEDIGREE ""
    set L_others ""
    set L_CHROM ""
    set F [open $tmpOUTPUTFILE]
    while {[gets $F L]>=0} {
        if {[regexp "^#" $L]} {
            if {[regexp "^##fileformat=" $L]} {
                lappend L_fileformat $L
                continue
            } elseif {[regexp "^##assembly=" $L]} {
                lappend L_assembly "$L (before using the liftoverSV)"
                continue
            } elseif {[regexp "^##contig=<ID" $L]} {
                lappend L_contig $L
                continue
            } elseif {[regexp "^##FILTER=<ID" $L]} {
                lappend L_FILTER $L
                continue
            } elseif {[regexp "^##INFO=<ID" $L]} {
                lappend L_INFO $L
                continue
            } elseif {[regexp "^##FORMAT=<ID" $L]} {
                lappend L_FORMAT $L
                continue
            } elseif {[regexp "^##ALT=<ID" $L]} {
                lappend L_ALT $L
                continue
            } elseif {[regexp "^##META=<ID" $L]} {
                lappend L_META $L
                continue
            } elseif {[regexp "^##SAMPLE=<ID" $L]} {
                lappend L_SAMPLE $L
                continue
            } elseif {[regexp "^##PEDIGREE=<ID" $L]} {
                lappend L_PEDIGREE $L
                continue
            } elseif {[regexp "^#CHROM" $L]} {
                set L_CHROM $L
            } else {
                lappend L_others $L
            }
        } else {
            break
        }
    }
    
    
    set L_toWrite {}
    
    ##fileformat=VCFv4.4
    if {$L_fileformat ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_fileformat]}
    
    ##assembly=url
    # Indicate the use of liftoverSV in the header lines
    if {$L_assembly eq ""} {
        lappend L_assembly "##assembly=liftoverSV used with $g_liftoverSV(CHAIN)"
    }
    lappend L_assembly "##liftoverSV_version=$g_liftoverSV(Version); $g_liftoverSV(CHAIN)"
    lappend L_toWrite {*}[lsort -dictionary $L_assembly]
    
    ##contig=<ID=ctg1,length=81195210,URL=ftp://somewhere.example/assembly.fa,md5=f126cdf8a6e0c7f379d618ff66beb2da,...>
    if {$L_contig ne ""} {
        foreach contigLine [lsort -dictionary $L_contig] {
            # Add the length of the contig if absent:
            if {![regexp "length=" $contigLine]} {
                if {[regexp "^##contig=<ID=(\[^,>\]+)" $contigLine match contig]} {
                    catch {regsub "^##contig=<ID=$contig" $contigLine "##contig=<ID=$contig, length=$g_liftoverSV(sizeAfterLift,$contig)" contigLine} Message
                }
            }
            lappend L_toWrite $contigLine
        }
    }
    
    if {$L_FILTER ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_FILTER]}
    if {$L_INFO ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_INFO]}
    if {$L_FORMAT ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_FORMAT]}
    if {$L_ALT ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_ALT]}
    if {$L_META ne ""} {lappend L_toWrite {*}[join [lsort -dictionary $L_META]}
    if {$L_SAMPLE ne ""} {lappend L_toWrite {*}[join [lsort -dictionary $L_SAMPLE]}
    if {$L_PEDIGREE ne ""} {lappend L_toWrite {*}[join [lsort -dictionary $L_PEDIGREE]}
    if {$L_others ne ""} {lappend L_toWrite {*}[lsort -dictionary $L_others]}
    lappend L_toWrite $L_CHROM
    
    
    WriteTextInFile [join $L_toWrite "\n"] $g_liftoverSV(OUTPUTFILE)
    eval exec grep -v "^#" $tmpOUTPUTFILE | sort -k1,1V -k2,2n >> $g_liftoverSV(OUTPUTFILE)
    
    
    file delete -force $tmpOUTPUTFILE
    
    
    return
}


proc sortTheLiftedVCF {} {
    
    global g_liftoverSV
    
    
    regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".sort.vcf.gz" sortedOUTPUTFILE
    puts "\n...writing [file tail $sortedOUTPUTFILE] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    
    # Try to sort the output file
    if {[catch {eval exec $g_liftoverSV(BCFTOOLS) sort $g_liftoverSV(OUTPUTFILE) -o $sortedOUTPUTFILE -Oz} Message]} {
        if {![regexp "Done$" $Message]} {
            puts "bcftools sort $g_liftoverSV(OUTPUTFILE) -o $sortedOUTPUTFILE -Oz"
            puts $Message
        }
    }
    
    # Clean
    if {[file exists $sortedOUTPUTFILE]} {
        file delete -force $g_liftoverSV(OUTPUTFILE)
    }
    
}

