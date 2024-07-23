############################################################################################################
# liftoverSV 0.1.1_beta                                                                                    #
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

proc configureLiftoverSV {argv} {

    global g_liftoverSV

    puts "...loading configuration data ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"


    #######################
    ## Load default options
    #######################
    set g_liftoverSV(BEDTOOLS)      "bedtools"
	set g_liftoverSV(LIFTOVER)      "liftOver"   
    set g_liftoverSV(PERCENT)       "0.05"


    ##################################
    ## Load options given in arguments
    ##################################
    set lOptionsOk "B BEDTOOLS C CHAIN h help I INPUTFILE L LIFTOVER O OUTPUTFILE P PERCENT R REFFASTASEQ"

    set i 0
    set j 1
    while {$j < [llength $argv]} {
        set optionName [lindex $argv $i]
        regsub "^-(-?)" $optionName "" optionName
        set optionValue [lindex $argv $j]
        set  k [lsearch -exact $lOptionsOk $optionName]
        if {$k != -1} {
			if {$optionName eq "B"} {
				set optionName BEDTOOLS
			} elseif {$optionName eq "C"} {
                set optionName CHAIN
            } elseif {$optionName eq "I"} {
                set optionName INPUTFILE
            } elseif {$optionName eq "L"} {
                set optionName LIFTOVER
            } elseif {$optionName eq "O"} {
                set optionName OUTPUTFILE
            } elseif {$optionName eq "P"} {
                set optionName PERCENT
            } elseif {$optionName eq "R"} {
                set optionName REFFASTASEQ
            } elseif {$optionName eq "h" || $optionName eq "help"} {
				showHelp; exit 0
			} 
            set g_liftoverSV($optionName) $optionValue
        } else {
            puts "\"$optionName\" option not known."
            puts "For more information on the arguments, please use the --help option"
            exit 2
        }
        incr i 2
        incr j 2
    }

    ########################################
    ## Checking of the configuration options
    ########################################
    puts "...checking configuration data ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

    ## INPUTFILE: We should have a VCF input file
    if {$g_liftoverSV(INPUTFILE) eq ""} {
        puts "LiftoverSV needs in argument the path of your SV VCF input file (--INPUTFILE ...) - Exit with error."
        exit 2
    }
    if {![regexp -nocase "\\.(vcf(.gz)?)$" $g_liftoverSV(INPUTFILE)]} {
        puts "############################################################################"
        puts "Bad option value: --INPUTFILE = $g_liftoverSV(INPUTFILE)"
        puts "Extension file should be \".vcf\" - Exit with error."
        puts "############################################################################"
        exit 2
    }
    ## INPUTFILE: It must be an existing file
    if {![file exists $g_liftoverSV(INPUTFILE)]} {
        puts "############################################################################"
        puts "Bad value for the --INPUTFILE option, file does not exist ($g_liftoverSV(INPUTFILE)) - Exit with error."
        puts "############################################################################"
        exit 2
    }
    ## INPUTFILE: no SV to lift if it is an empty file
    if {[isAnEmptyFile $g_liftoverSV(INPUTFILE)]} {
        puts "############################################################################"
        puts "INPUTFILE ($g_liftoverSV(INPUTFILE) is empty, no SV to lift - Exit without error."
        puts "############################################################################"
        exit 0
    }

	## CHAIN: We should have a ".chain" file
    if {$g_liftoverSV(CHAIN) eq ""} {
        puts "LiftoverSV needs in argument the path of your CHAIN file (--CHAIN ...) - Exit with error."
        exit 2
    }
    if {![regexp -nocase "(\\.chain)$" $g_liftoverSV(CHAIN)]} {
        puts "############################################################################"
        puts "Bad option value: --CHAIN = $g_liftoverSV(CHAIN)"
        puts "Extension file should be \".chain\" - Exit with error."
        puts "############################################################################"
        exit 2
    }
    ## CHAIN: It must be an existing file
    if {![file exists $g_liftoverSV(CHAIN)]} {
        puts "############################################################################"
        puts "Bad value for the --CHAIN option, file does not exist ($g_liftoverSV(CHAIN)) - Exit with error."
        puts "############################################################################"
        exit 2
    }

    ## REFFASTASEQ: If defined, we should have an existing ".fasta" file
    if {![info exists g_liftoverSV(REFFASTASEQ)]} {
        #puts "LiftoverSV needs in argument the path of your REFFASTASEQ file (--REFFASTASEQ ...) - Exit with error."
        #exit 2
    } elseif {![regexp -nocase "(\\.fasta)$" $g_liftoverSV(REFFASTASEQ)]} {
        puts "############################################################################"
        puts "Bad option value: --REFFASTASEQ = $g_liftoverSV(REFFASTASEQ)"
        puts "Extension file should be \".fasta\" - Exit with error."
        puts "############################################################################"
        exit 2
    } elseif {![file exists $g_liftoverSV(REFFASTASEQ)]} {
        puts "############################################################################"
        puts "Bad value for the --REFFASTASEQ option, file does not exist ($g_liftoverSV(REFFASTASEQ)) - Exit with error."
        puts "############################################################################"
        exit 2
    }

    ## OUTPUTFILE extension must be ".vcf".
    if {$g_liftoverSV(OUTPUTFILE) eq ""} {
        puts "LiftoverSV needs in argument the path of your OUTPUTFILE file (--OUTPUTFILE ...) - Exit with error."
        exit 2
    }
    if {![regexp -nocase "(\\.vcf)$" $g_liftoverSV(OUTPUTFILE)]} {
        puts "############################################################################"
        puts "Bad option value: --OUTPUTFILE = $g_liftoverSV(OUTPUTFILE)"
        puts "Extension file should be \".vcf\" - Exit with error."
        puts "############################################################################"
        exit 2
    }

    # BEDTOOLS, LIFTOVER: It should be a good path that we can run
    foreach tool {BEDTOOLS LIFTOVER} {
        if {[catch {eval exec $g_liftoverSV($tool)} Message] && ![regexp "usage:" $Message]} {
            puts "############################################################################"
            puts "Bad value for the $tool option ($g_liftoverSV($tool))"
            puts "$Message"
            puts "Exit with error."
            puts "############################################################################"
            exit 2
        } 
    }

    ## PERCENT: should be defined between 0 and 1.
    if {$g_liftoverSV(PERCENT) < 0 || $g_liftoverSV(PERCENT) > 1} {
        puts "Bad option value: --PERCENT = $g_liftoverSV(PERCENT)"
        puts "Should be in the \[0-1\] range values, default = 0.05"
        exit 2
    }

 
}
