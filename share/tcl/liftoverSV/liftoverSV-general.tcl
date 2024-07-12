############################################################################################################
# LiftOver_SVs_from_VCF 1.0.0_beta                                                                         #
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

