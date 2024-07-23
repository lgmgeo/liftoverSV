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

proc showHelp {} {

    global g_liftoverSV

    if {[file exists $g_liftoverSV(installDir)/commandLineOptions.txt]} {
        puts [ContentFromFile $g_liftoverSV(installDir)/commandLineOptions.txt]
    } elseif {[file exists $g_liftoverSV(docDir)/commandLineOptions.txt]} {
        puts [ContentFromFile $g_liftoverSV(docDir)/commandLineOptions.txt]
    }

    return
}

