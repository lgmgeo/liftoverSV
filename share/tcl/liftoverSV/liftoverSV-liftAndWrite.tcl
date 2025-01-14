############################################################################################################
# liftoverSV 0.1.2_beta                                                                                    #
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


# Extract the different genomic coordinates from the input VCF
proc extractAllTheGenomicCoordinates {} {

	global g_liftoverSV


	puts "\n...reading [file tail $g_liftoverSV(INPUTFILE)] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	set tmpBEDfile "./[clock format [clock seconds] -format "%Y%m%d-%H%M%S"].tmp.bed"
	puts "\t* creation of $tmpBEDfile"

	set L_toWrite {}
	set i 0

	if {[regexp ".gz$" $g_liftoverSV(INPUTFILE)]} {
	    set F [open "|gzip -cd $g_liftoverSV(INPUTFILE)"]
	} else {
	    set F [open "$g_liftoverSV(INPUTFILE)"]
	}
	while {[gets $F L]>=0} {
	
		if {[regexp "^#" $L]} {continue}
		set Ls [split $L "\t"]

	    set chrom [lindex $Ls 0]
	    set start [lindex $Ls 1]
		set alt [lindex $Ls 4]
		set infos [lindex $Ls 7]

	    lappend L_toWrite "$chrom\t$start\t[expr {$start+1}]\t$i"
		incr i

		# INFO/END
	    if {[regexp "(^END|;END)=(\[^;\]+)(;|$)" $infos match tutu end titi]} {
			if {![info exists localMemory($chrom,$end)]} {
				set localMemory($chrom,$end) 1
				lappend L_toWrite "$chrom\t$end\t[expr {$end+1}]\t$i"
				incr i
			}
		}

		# INFO/SVEND
		if {[regexp "(^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos match tutu svend titi]} {
            if {![info exists localMemory($chrom,$svend)]} {
                set localMemory($chrom,$svend) 1
				lappend L_toWrite "$chrom\t$svend\t[expr {$svend+1}]\t$i"
				incr i
			}
		}

		# Square-bracketed ALT notation
		if {[regexp "(\\\[|\\\])(\[chr0-9\]+):(\[0-9\]+)(\\\[|\\\])" $alt match tutu chromALT posALT]} {
			if {![info exists localMemory($chromALT,$posALT)]} {
                set localMemory($chromALT,$posALT) 1
				lappend L_toWrite "$chromALT\t$posALT\t[expr {$posALT+1}]\t$i"
				incr i
			}
        }


		if {![expr {$i%100000}]} {
			WriteTextInFile [join $L_toWrite "\n"] $tmpBEDfile
			set L_toWrite {}
		}
	}
	close $F
	WriteTextInFile [join $L_toWrite "\n"] $tmpBEDfile
	
	return $tmpBEDfile
}



# Liftover these different genomic coordinates
proc liftoverGenomicCoordinates {tmpBEDfile} {

    global g_liftoverSV

	puts "\n...liftover $tmpBEDfile ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	if {![catch {set i [lindex [eval exec wc -l $tmpBEDfile] 0]} Message]} {
		puts "\t($i genomic coordinates)"
	}

	regsub ".bed" $tmpBEDfile ".lifted.bed" tmpBEDfileLifted
	set command "$g_liftoverSV(LIFTOVER) $tmpBEDfile $g_liftoverSV(CHAIN) $tmpBEDfileLifted unMapped"
	puts "\t* Command line:"
	puts "\t\t$command"

	if {[catch {eval exec $command} Message]} {
		puts "\t* Output:"
		foreach L [split $Message "\n"] {
			puts "\t\t$L"
		}
	}

	file delete -force unMapped
	return $tmpBEDfileLifted
}



# Memorize the lifted coordinates
proc memorizeGenomicCoordinates {tmpBEDfile tmpBEDfileLifted} {

    global g_liftoverSV
	global g_ID
	global g_liftedCoord
	global g_L_chrom

	puts "\n...memorizing the lifted coordinates ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    set g_L_chrom(notLifted) {}
	set F [open "$tmpBEDfile"]
	while {[gets $F L]>=0} {
	    set Ls [split $L "\t"]
		set chrom [lindex $Ls 0]
		set pos [lindex $Ls 1]
	    set id [lindex $Ls end]
	    set g_ID($chrom,$pos) $id
        if {[lsearch -exac $g_L_chrom(notLifted) $chrom] eq -1} {
            lappend g_L_chrom(notLifted) $chrom
        }
	}
	file delete -force $tmpBEDfile

	set g_L_chrom(lifted) {}
	set F [open "$tmpBEDfileLifted"]
	while {[gets $F L]>=0} {
	    set Ls [split $L "\t"]
		set id [lindex $Ls end]
        set chrom [lindex $Ls 0]
        set pos [lindex $Ls 1]
		set g_liftedCoord($id) "$chrom,$pos"
		if {[lsearch -exac $g_L_chrom(lifted) $chrom] eq -1} {
			lappend g_L_chrom(lifted) $chrom
		}
	}
	file delete -force $tmpBEDfile
	file delete -force $tmpBEDfileLifted
}



# Write the lifted VCF
# Rules:
# => lift the CHROM, the END/SVEND and the REF/ALT of the SV.
# => drop the SV if:
#   - Case1: one position (start or end) is lifted while the other doesn't
#   - Case2: one position (start or end) goes to a different chrom from the other
#   - Case3: "lifted start" > "lifted end"
#   - Case4: the distance between the two lifted positions changes significantly (difference between both SVLENs > 5%)

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
# INFO/CIPOS and INFO/CIEND:
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
        puts "Exit with error."
        exit 2
    }

	puts "\n...writing [file tail $g_liftoverSV(OUTPUTFILE)] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	file delete -force $g_liftoverSV(OUTPUTFILE)
    file delete -force $unmappedFile

	# SVLEN calculation only for the following svtypes (not INS, TRA, CPX...)
	set L_svtypes {DUP DEL INV}

	# Memorize if there are new INFO/FORMAT annotations to add in the header lines
    set L_header_INFO ""
    set L_header_FORMAT ""
	set L_SVlines_INFO ""
	set L_SVlines_FORMAT ""
    set L_new_INFO ""
    set L_new_FORMAT ""

	set L_toWrite {}
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
	set testContigs 1
	while {[gets $F L]>=0} {

		# VCF header lines
		##################
	    if {[regexp "^#" $L]} {
			if {[regexp "^##contig=<ID=" $L]} {
				##contig=<ID=chr1,length=249250621>
				regsub ",length=\[0-9\]+" $L "" L
				set testContigs 2
				lappend L_toWrite "$L"
				continue
			}
			if {$testContigs eq 2} {
				foreach chrom [lsort $g_L_chrom(lifted)] {
					if {[lsearch -exact $g_L_chrom(notLifted) $chrom] ne -1} {
						lappend L_toWrite "##contig=<ID=$chrom>"
					}
				}
				set testContigs 0
			}
            lappend L_toWrite "$L"

			##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">
			##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
			if {[regexp "^##INFO=<ID=(\[^,\]+)," $L match knownValue]} {
				lappend L_header_INFO $knownValue
			}
            if {[regexp "^##FORMAT=<ID=(\[^,\]+)," $L match knownValue]} {
                lappend L_header_FORMAT $knownValue
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

        set svtype [normalizeSVtype $chrom $start $ref $alt]

		# Drop SV if the SVTYPE is not DUP, DEL, INV, INS or TRA
		# if {$svtype eq "None"} {continue}
		# ...A voir si nécessaire

		# Drop SV without a square or an angle bracketed notation.
        # ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
		# Case5
        if {[regexp -nocase "^\[acgtnACGTN.*\]+$" $ref$alt]} {
            lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tnot a square or angle bracketed notation for the ALT feature."
			incr n_unmapped 
			incr case5
            continue
        }

		# Memorize all the INFO annotations presents in the SV lines
		regsub -all "=\[^;\]+" $infos "" infosToCheck
		lappend L_SVlines_INFO {*}[split $infosToCheck ";"]
		set L_SVlines_INFO [lsort -unique $L_SVlines_INFO]

        # Memorize all the FORMAT annotations presents in the SV lines
        set format_annot [lindex $Ls 8]
        lappend L_SVlines_FORMAT {*}[split $format_annot ";"]
        set L_SVlines_FORMAT [lsort -unique $L_SVlines_FORMAT]


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
	

		# Lift of the sequence in REF
        if {[regexp "^N+$" $ref]} {
			set theNewRef $ref       ;# REF is composed only of "N"
		} elseif {$ref eq "."} {
            set theNewRef "."        ;# REF = "."
		} else {
            regsub -all "\[*.\]" $ref "" refbis
	        set REFlength [string length $refbis] ;# REF = "A" or "CTTGA" or ...
			set theNewRef [ExtractDNAseq $theNewChrom [expr {$theNewStart-1}] [expr {$theNewStart+$REFlength-1}]]
		}


		# Lift of the square-bracketed ALT 
		set theNewALT $alt
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
			# A]chr2:32156] => "A" = $theNewRef
			# [chr2:32156[A => "A" = $theNewRef
			# INS:
			# ACCCCC[chr2:32156[ => chr2:32156 is "A" (first NT before the inserted sequence)
			# ]chr2:32156]CCCCCA => chr2:32156 is "A"(last NT after the inserted sequence)
			if {$baseLeft ne "" && $baseRight eq ""} {
				set theNewALT "$theNewRef[string range $baseLeft 1 end]$bracketLeft$theNewChromALT:$theNewPosALT$bracketRight"
			} elseif {$baseRight ne "" && $baseLeft eq ""} {
				set theNewALT "$bracketLeft$theNewChromALT:$theNewPosALT$bracketRight[string range $baseRight 0 end-1]$theNewRef"
			} else {
				lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe ALT features ($alt) is badly formatted."
				incr n_unmapped
				incr case5
				continue
			}
		}


		# Memorize the 6 first columns lifted
		set theNewL "$theNewChrom\t$theNewStart\t$ID\t$theNewRef\t$theNewALT\t$qual\t$filter"


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
		# Check and modify if needed the CIPOS/CIEND values in order to have:
		#   => POS-CIPOS > 0
		#   => END+CIEND < chrom_length
		if {[regexp "(^CIPOS|;CIPOS)=(\[^;\]+)(;|$)" $infos match titi cipos]} {
			if {[expr {$theNewStart-$cipos}] <= 0} {
				set cipos [expr {$theNewStart-1}]
				regsub "(^CIPOS|;CIPOS)=(\[0-9\]+)(;|$)" $infos "\\1=$cipos\\3" infos
			}
		}
        if {[regexp "(^CIEND|;CIEND)=(\[^;\]+)(;|$)" $infos match titi ciend]} {
            if {[expr {$theNewEND+$ciend}] > $g_liftoverSV(sizeAfterLift,$theNewChrom)} {
				set ciend [expr {$g_liftoverSV(sizeAfterLift,$theNewChrom)-$theNewEND}]
				regsub "(^CIEND|;CIEND)=(\[0-9\]+)(;|$)" $infos "\\1=$ciend\\3" infos
			}
        }


		# tmp output (before the sort)
		##############################
		append theNewL "\t$infos\t[join [lrange $Ls 8 end] "\t"]"

		lappend L_toWrite $theNewL
		incr n_mapped

	    if {![expr {$n_mapped%10000}]} {
	        WriteTextInFile [join $L_toWrite "\n"] $g_liftoverSV(OUTPUTFILE)
	        set L_toWrite {}
	    }
        if {$n_unmapped ne 0 && ![expr {$n_unmapped%10000}]} {
            WriteTextInFile [join $L_unmappedToWrite "\n"] $unmappedFile
            set L_unmappedToWrite {}
        }

	}
	close $F

	# Writings
	if {$L_toWrite ne ""} {
		WriteTextInFile [join $L_toWrite "\n"] $g_liftoverSV(OUTPUTFILE)
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

	# Add new header lines (only for INFO or FORMAT) if needed
	foreach val $L_SVlines_INFO {
		if {[lsearch -exact $L_header_INFO $val] eq -1} {lappend L_new_INFO $val}
	}
    foreach val $L_SVlines_FORMAT {
        if {[lsearch -exact $L_header_FORMAT $val] eq -1} {lappend L_new_FORMAT $val}
    }

	if {$L_new_INFO ne "" || $L_new_FORMAT ne ""} {
		addNewHeaderLinesAndSort $L_new_INFO $L_new_FORMAT
	}

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

# Add new header lines (only for INFO or FORMAT) if needed:
###########################################################
# As the format (Number, String) is not known, "Number=." and "Type=String" values are used by default:
# - ##FORMAT=<ID=XXX,Number=.,Type=String,Description="XXX">
# - ##INFO=<ID=YYY,Number=.,Type=String,Description="YYY">
proc addNewHeaderLinesAndSort {L_new_INFO L_new_FORMAT} {

    global g_liftoverSV
 
    regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".checked.vcf" checkedOUTPUTFILE
	file delete -force $checkedOUTPUTFILE

	foreach val $L_new_INFO {
		lappend L_INFO "##INFO=<ID=$val,Number=.,Type=String,Description=\"$val\">"
	}
    foreach val $L_new_FORMAT {
        lappend L_FORMAT "##FORMAT=<ID=$val,Number=.,Type=String,Description=\"$val\">"
    }

    set L_fileformat ""
    set L_assembly ""
	set L_contig ""
	set L_FILTER ""
	set L_ALT ""
	set L_META ""
	set L_SAMPLE ""
	set L_PEDIGREE ""
	set L_others ""
	set L_CHROM ""
	set F [open $g_liftoverSV(OUTPUTFILE)]
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
	
	if {$L_fileformat ne ""} {WriteTextInFile [join [lsort -dictionary $L_fileformat] "\n"] $checkedOUTPUTFILE}

	# Indicate the use of liftoverSV in the header lines
    if {$L_assembly eq ""} {lappend L_assembly "##assembly=liftoverSV used with $g_liftoverSV(CHAIN)"}
	WriteTextInFile [join [lsort -dictionary $L_assembly] "\n"] $checkedOUTPUTFILE
    set newHeaderLine "##liftoverSV_version=$g_liftoverSV(Version); $g_liftoverSV(CHAIN); [clock format [clock seconds] -format "%B %d %Y %H:%M"]"
    WriteTextInFile $newHeaderLine $checkedOUTPUTFILE

    if {$L_contig ne ""} {WriteTextInFile [join [lsort -dictionary $L_contig] "\n"] $checkedOUTPUTFILE}
    if {$L_FILTER ne ""} {WriteTextInFile [join [lsort -dictionary $L_FILTER] "\n"] $checkedOUTPUTFILE}
    if {$L_INFO ne ""} {WriteTextInFile [join [lsort -dictionary $L_INFO] "\n"] $checkedOUTPUTFILE}
    if {$L_FORMAT ne ""} {WriteTextInFile [join [lsort -dictionary $L_FORMAT] "\n"] $checkedOUTPUTFILE}
    if {$L_ALT ne ""} {WriteTextInFile [join [lsort -dictionary $L_ALT] "\n"] $checkedOUTPUTFILE}
    if {$L_META ne ""} {WriteTextInFile [join [lsort -dictionary $L_META] "\n"] $checkedOUTPUTFILE}
    if {$L_SAMPLE ne ""} {WriteTextInFile [join [lsort -dictionary $L_SAMPLE] "\n"] $checkedOUTPUTFILE}
    if {$L_PEDIGREE ne ""} {WriteTextInFile [join [lsort -dictionary $L_PEDIGREE] "\n"] $checkedOUTPUTFILE}
    if {$L_others ne ""} {WriteTextInFile [join [lsort -dictionary $L_others] "\n"] $checkedOUTPUTFILE}

	WriteTextInFile $L_CHROM $checkedOUTPUTFILE

	eval exec grep -v "^#" $g_liftoverSV(OUTPUTFILE) | sort -k1,1V -k2,2n >> $checkedOUTPUTFILE

	file rename -force $checkedOUTPUTFILE $g_liftoverSV(OUTPUTFILE)

	return
}


proc sortTheLiftedVCF {} {

    global g_liftoverSV

	# Try to sort the output file
	regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".sort.vcf.gz" sortedOUTPUTFILE

	if {[catch {eval exec $g_liftoverSV(BCFTOOLS) sort $g_liftoverSV(OUTPUTFILE) -o $sortedOUTPUTFILE -Oz} Message]} {
		if {![regexp "Done$" $Message]} {
			puts "\n\nPlease sort and compress the output file, run :"
			puts "bcftools sort $g_liftoverSV(OUTPUTFILE) -o $sortedOUTPUTFILE -Oz"
			puts $Message
		} else {
			puts "\n...sort and compress [file tail $g_liftoverSV(OUTPUTFILE)] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
			puts "\t* creation of [file tail $sortedOUTPUTFILE]"
			if {[file exists $sortedOUTPUTFILE]} {
				file delete -force $g_liftoverSV(OUTPUTFILE)
			}
		}
	}
}


