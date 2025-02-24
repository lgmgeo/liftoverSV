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


# Read a .chain file to retrieve the size of the chromosomes after the liftover
proc retrieveChromSize {} {
    
    global g_liftoverSV
    
    puts "\n...memorizing the size of the chromosomes for the query build ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
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



# Check if the chain file and the reference fasta file are coherent:
# - both with or without "chr"
proc checkREFFASTASEQ {} {
    
    global g_liftoverSV
    
    # chain 20851231461 chr1 249250621 + 10000 249240621 chr1 248956422 + 10000 248946422 2
    
    puts "\n...checking the REFFASTASEQ file ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    # Check the CHAIN file
    set f [open "$g_liftoverSV(CHAIN)"]
    while {![eof $f]} {
        set L [gets $f]
        if {[regexp "^chain" $L]} {
            if {[regexp "chr" [lindex $L 7]]} {set CHAINstatus "with"} else {set CHAINstatus "without"}
            break
        }
    }
    close $f
    
    # Check the REFFASTASEQ
    set f [open $g_liftoverSV(REFFASTASEQ)]
    while {![eof $f]} {
        set L [gets $f]
        if {[regexp "^>" $L]} {
            if {[regexp "chr" [lindex $L 0]]} {set FASTAstatus "with"} else {set FASTAstatus "without"}
            break
        }
    }
    close $f
    
    # Compare both files
    if {$CHAINstatus ne $FASTAstatus} {
        puts "\nIncoherence:"
        puts "############"
        puts "- CHAIN: $g_liftoverSV(CHAIN)"
        puts "     => Contig names of the query build: $CHAINstatus the 'chr' prefix"
        puts "- REFFASTASEQ: $g_liftoverSV(REFFASTASEQ)"
        puts "     => Contig names of the query build: $FASTAstatus the 'chr' prefix"
        puts "\nPlease, check your REFFASTASEQ option"
        puts "\nExit"
        exit 2
    }
    
    return
}


# Extract DNA sequences from genomic coordinates
################################################
# Input = 1-based coordinates (VCF)
#
# Equivalence between BED and VCF coordinates:
#          VCF [start, end]
#          BED [start-1, end)
#
# Return:
# - the sequence
# - "" in case of error
#
# Example: VCFchrom=chr22 VCFstart=16848507 VCFend=16848507
#		   => g (hg38)
proc ExtractDNAseq {VCFchrom VCFstart VCFend} {
    
    global g_liftoverSV
    
    # Use of "bedtools getfasta"
    ############################
    # Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>
    # Options:
    #          -fi			Input FASTA file
    #          -bed			BED/GFF/VCF file of ranges to extract from -fi
    #          -fo			Output file (can be FASTA or TAB-delimited)
    
    set outputSeq ""
    set BEDstart [expr {$VCFstart-1}]
    set command "echo -e \"$VCFchrom\t$BEDstart\t$VCFend\" | bedtools getfasta -fi $g_liftoverSV(REFFASTASEQ) -bed -"
    if {[catch {eval exec $command} Message]} {
        #puts $command
        #puts $Message
    } else {
        foreach L [split $Message "\n"] {
            if {[regexp "^>" $L]} {continue}
            if {$L eq ""} {continue}
            append outputSeq $L
        }
    }
    
    return $outputSeq
}



# Check if the VCF input file contains multi-allelic lines
proc isMultiAllelic {} {
    
    global g_liftoverSV
    
    
    if {[regexp -nocase ".gz$" $g_liftoverSV(INPUTFILE)]} {
        set command "zcat $g_liftoverSV(INPUTFILE) | grep -v ^# | cut -f 4-5  | grep -c ,"
    } else {
        set command "grep -v ^# $g_liftoverSV(INPUTFILE) | cut -f 4-5  | grep -c ,"
    }
    
    catch {eval exec $command} nMultiAllelicLines
    
    if {[string index $nMultiAllelicLines 0] ne "0"} {
        puts "####################################################################################"
        puts "Please split the multi-allelic lines of the VCF input file before to run liftoverSV"
        puts "Exit without error."
        puts "####################################################################################"
        exit 0
    }
    
    return
}

