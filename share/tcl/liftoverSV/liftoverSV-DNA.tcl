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


# Read a .chain file to retrieve the size of the chromosomes before / after the liftover
proc retrieveChromSize {} {

    global g_liftoverSV

	puts "...memorizing the size of the chromosomes for the query build"
	foreach L [LinesFromFile $g_liftoverSV(CHAIN)] {
		if {[regexp "^chain" $L]} {
			#set chromBeforeLift [lindex $L 2]
			#set sizeBeforeLift [lindex $L 3]
            #set g_liftoverSV(sizeBeforeLift,$chromBeforeLift) $sizeBeforeLift
            set chromAfterLift [lindex $L 7]
            set sizeAfterLift [lindex $L 8]
            set g_liftoverSV(sizeAfterLift,$chromAfterLift) $sizeAfterLift
		}
	}
	return
}


# Extract DNA sequences into a fasta file based on feature coordinates
proc ExtractDNAseq {chrom start end FastaFile} {
 
	global g_liftoverSV

	# Use of "bedtools getfasta"
	############################
	# Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>
	# Options:
	#          -fi			Input FASTA file
	#          -bed			BED/GFF/VCF file of ranges to extract from -fi
	#          -fo			Output file (can be FASTA or TAB-delimited)

	set tmpCoordFile "./[clock format [clock seconds] -format "%Y%m%d-%H%M%S"].tmp.coord.bed"
	WriteTextInFile "$chrom\t$start\t$end" $tmpCoordFile 

	set tmpOutputFile "./[clock format [clock seconds] -format "%Y%m%d-%H%M%S"].tmp.coord.fasta"

	set command "$g_liftoverSV(bedtools) getfasta -fi $FastaFile -bed $tmpCoordFile -fo $tmpOutputFile"
	if {[catch {eval exec $command} Message]} {
		puts $command
		puts $Message
	}

    set f [open $tmpOutputFile]
    while {![eof $f]} {
        set L [gets $f]
		if {[regexp "^>" $L]} {continue}
		if {$L eq ""} {continue}
		append outputSeq $L
	}
	close $f

	file delete -force $tmpCoordFile
	file delete -force $tmpOutputFile

	return $outputSeq
}





