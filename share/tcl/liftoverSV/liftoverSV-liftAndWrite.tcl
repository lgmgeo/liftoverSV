############################################################################################################
# liftoverSV 1.0.0_beta                                                                                    #
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
	puts "   => creation of $tmpBEDfile"

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
		#set alt [lindex $Ls 4]
		set infos [lindex $Ls 7]

	    lappend L_toWrite "$chrom\t$start\t[expr {$start+1}]\t$i"
		incr i

	    if {[regexp "(^END|;END)=(\[^;\]+)(;|$)" $infos match tutu end titi]} {
		    lappend L_toWrite "$chrom\t$end\t[expr {$end+1}]\t$i"
			incr i
		}

		if {[regexp "(^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos match tutu svend titi]} {
			lappend L_toWrite "$chrom\t$svend\t[expr {$svend+1}]\t$i"
			incr i
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
	regsub ".bed" $tmpBEDfile ".lifted.bed" tmpBEDfileLifted
	set command "$g_liftoverSV(LIFTOVER) $tmpBEDfile $g_liftoverSV(CHAIN) $tmpBEDfileLifted unMapped"
	puts $command

	if {[catch {eval exec $command} Message]} {
		puts $Message
	}

	file delete -force unMapped
	return $tmpBEDfileLifted
}



# Memorize the lifted coordinates
proc memorizeGenomicCoordinates {tmpBEDfile tmpBEDfileLifted} {

    global g_liftoverSV
	global g_ID
	global g_liftedCoord

	puts "\n...memorizing the lifted coordinates ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	set F [open "$tmpBEDfile"]
	while {[gets $F L]>=0} {
	    set Ls [split $L "\t"]
	    set id [lindex $Ls end]
	    set g_ID([lindex $Ls 0],[lindex $Ls 1]) $id
	}
	file delete -force $tmpBEDfile

	set F [open "$tmpBEDfileLifted"]
	while {[gets $F L]>=0} {
	    set Ls [split $L "\t"]
		set id [lindex $Ls end]
		set g_liftedCoord($id) "[lindex $Ls 0],[lindex $Ls 1]"
	}
	file delete -force $tmpBEDfile
	file delete -force $tmpBEDfileLifted
}



# Write the lifted VCF
# Rules:
# => lift the CHROM and the END/SVEND of the SV.
# => drop the SV if:
#   - Case1: one position (start or end) is lifted while the other doesn't
#   - Case2: one position (start or end) goes to a different chrom from the other
#   - Case3: "lifted start" > "lifted end"
#   - Case4: the distance between the two lifted positions changes significantly (difference between both SVLENs > 5%)
proc writeTheLiftedVCF {} {

    global g_liftoverSV
    global g_ID
    global g_liftedCoord


    if {![regsub ".vcf$" $g_liftoverSV(OUTPUTFILE) ".unmapped" unmappedFile]} {
        puts "Check the output file name ($g_liftoverSV(OUTPUTFILE))."
        puts "Extension sould be \".vcf\""
        puts "Exit with error."
        exit 2
    }

	puts "\n...writing [file tail $g_liftoverSV(OUTPUTFILE)] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	file delete -force $g_liftoverSV(OUTPUTFILE)
    file delete -force $unmappedFile

	set L_toWrite {}
	set i 1
	set L_unmappedToWrite {}
	set j 1

	set case1 0
    set case2 0
    set case3 0
    set case4 0

	if {[regexp ".gz$" $g_liftoverSV(INPUTFILE)]} {
	    set F [open "|gzip -cd $g_liftoverSV(INPUTFILE)"]
	} else {
		set F [open "$g_liftoverSV(INPUTFILE)"]
	}
	while {[gets $F L]>=0} {

	    if {[regexp "^#" $L]} {
			regsub ",length=\[0-9\]+" $L "" L
			lappend L_toWrite "$L"
			continue
		}

		set Ls [split $L "\t"]

	    set chrom [lindex $Ls 0]
		set start [lindex $Ls 1]
		set infos [lindex $Ls 7]

		set theID $g_ID($chrom,$start)
		if {[info exists g_liftedCoord($theID)]} {
			set theNewChrom [lindex [split $g_liftedCoord($theID) ","] 0]
			set theNewStart [lindex [split $g_liftedCoord($theID) ","] 1]
		} else {
			# Case1
		    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tStart ($chrom,$start) not lifted"
			incr j
			incr case1
		    continue
		}
		set theNewL "$theNewChrom\t$theNewStart\t[join [lrange $Ls 2 6] "\t"]"

	    if {[regexp "(^END|;END)=(\[^;\]+)(;|$)" $infos match tutu end titi]} {
			set theID $g_ID($chrom,$end)
			if {[info exists g_liftedCoord($theID)]} {
	            set theNewChromEND [lindex [split $g_liftedCoord($theID) ","] 0]
				set theNewEND [lindex [split $g_liftedCoord($theID) ","] 1]
				if {$theNewChrom ne $theNewChromEND} {
					# Case2
					lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_start and lifted_END are located on different chrom ($theNewChrom # $theNewChromEND)"
					incr j
					incr case2
					continue
				}
				if {$theNewEND < $theNewStart} {
					# Case3
					lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_start ($theNewStart) > lifted_end ($theNewEND)"
					incr j
					incr case3
					continue
				}
				regsub "(^END|;END)=(\[^;\]+)(;|$)" $infos "\\1=$theNewEND\\3" infos
				set svlen [expr {$end-$start}]
                set svlenlifted [expr {$theNewEND-$theNewStart}]
				if {$svlenlifted < [expr {$svlen*0.95}] || $svlenlifted > [expr {$svlen*1.05}]} {
					# Case4
					lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe distance between the two lifted positions changes significantly (svlen diff > 5%)"
					incr j
                    incr case4
					continue
				}
			} else {
				# Case1
	            lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tEND ($chrom,$end) not lifted"
		        incr j
                incr case1
				continue
			}
	    }

	    if {[regexp "(^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos match tutu svend titi]} {
	        set theID $g_ID($chrom,$svend)
	        if {[info exists g_liftedCoord($theID)]} {
	            set theNewChromSVEND [lindex [split $g_liftedCoord($theID) ","] 0]
	            set theNewSVEND [lindex [split $g_liftedCoord($theID) ","] 1]
	            if {$theNewChrom ne $theNewChromSVEND} {
                    # Case2
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tstart and SVEND are located on different chrom ($theNewChrom # $theNewChromSVEND)"
                    incr j
                    incr case2
					continue
                }
				if {$theNewSVEND < $theNewStart} {
                    # Case3
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tlifted_start ($theNewStart) > lifted_svend ($theNewSVEND)"
                    incr j
                    incr case3
					continue
                }
	            regsub "(^SVEND|;SVEND)=(\[^;\]+)(;|$)" $infos "\\1=$theNewSVEND\\3" infos
                set svlen [expr {$svend-$start}]
				set svlenlifted [expr {$theNewSVEND-$theNewStart}]
                if {$svlenlifted < [expr {$svlen*0.95}] || $svlenlifted > [expr {$svlen*1.05}]} {
                    # Case4
                    lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tthe distance between the two lifted positions changes significantly (svlen diff > 5%; $svlen # $svlenlifted)"
                    incr j
                    incr case4
					continue
                }
            } else {
                # Case1
                lappend L_unmappedToWrite "[join [lrange $Ls 0 7] "\t"]\tSVEND ($chrom,$svend) not lifted"
                incr j
                incr case1
				continue
            }
	    }
	
	    if {[info exists svlenlifted]} {
			regsub "(^SVLEN|;SVLEN)=(-)?(\[0-9\]+)(;|$)" $infos "\\1=\\2$svlenlifted\\4" infos
	        regsub "(^SVSIZE|;SVSIZE)=(-)?(\[0-9\]+)(;|$)" $infos "\\1=\\2$svlenlifted\\4" infos
		} 
		catch {unset svlenlifted}

		append theNewL "\t$infos\t[join [lrange $Ls 8 end] "\t"]"

		lappend L_toWrite $theNewL
		incr i

	    if {![expr {$i%10000}]} {
	        WriteTextInFile [join $L_toWrite "\n"] $g_liftoverSV(OUTPUTFILE)
	        set L_toWrite {}
	    }
        if {![expr {$j%10000}]} {
            WriteTextInFile [join $L_unmappedToWrite "\n"] $unmappedFile
            set L_unmappedToWrite {}
        }

	}
	close $F
	WriteTextInFile [join $L_toWrite "\n"] $g_liftoverSV(OUTPUTFILE)
    WriteTextInFile [join $L_unmappedToWrite "\n"] $unmappedFile
	
	puts "\n...unmapped SV (see [file tail $unmappedFile] for details):"
	puts "\t- $case1 SV with one position (start or end) lifted while the other doesn't"
	puts "\t- $case2 SV with one position (start or end) mapped to a different chrom from the other"
	puts "\t- $case3 SV with \"lifted start\" > \"lifted end\""
	puts "\t- $case4 SV with the distance between the two lifted positions that changes significantly (difference between both SVLENs > 5%)"

	
}



