#!/usr/bin/env tclsh


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
# hg38 example: 
# ./extractDNAseq.tcl chr22 16848507 16848507 $REFFASTASEQ
# => g


set VCFchrom [lindex $argv 0]
set VCFstart [lindex $argv 1]
set VCFend [lindex $argv 2]
set REFFASTASEQ [lindex $argv 3]

# Use of "bedtools getfasta"
############################
# Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf> -fo <fasta>
# Options:
#          -fi          Input FASTA file
#          -bed         BED/GFF/VCF file of ranges to extract from -fi
#          -fo          Output file (can be FASTA or TAB-delimited)

set outputSeq ""
set BEDstart [expr {$VCFstart-1}]
set command "echo -e \"$VCFchrom\t$BEDstart\t$VCFend\" | bedtools getfasta -fi $REFFASTASEQ -bed -"
if {[catch {eval exec $command} Message]} {
    puts $command
    puts $Message
} else {
    foreach L [split $Message "\n"] {
        if {[regexp "^>" $L]} {continue}
        if {$L eq ""} {continue}
        append outputSeq $L
    }
}

puts $outputSeq


