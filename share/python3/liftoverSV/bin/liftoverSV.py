#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
liftoverSV 0.3.0_beta 
=====================

Copyright (C) 2024-current Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
import re
import time
import glob
import platform


def get_script_directory():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


def main(argv):

    # Initialisation of the g_liftoverSV dictionary
    ###############################################
    g_liftoverSV = {}

    # FHS directories setup:
    # - g_liftoverSV["install_dir"]
    # - g_liftoverSV["etc_dir"]
    # - g_liftoverSV["doc_dir"]
    # - g_liftoverSV["python_dir"]
    # - g_liftoverSV["bash_dir"]
    ######################################
    g_liftoverSV["install_dir"] = get_script_directory()
    # $LIFTOVERSV/share/python3/liftoverSV/bin/liftoverSV.py
    if g_liftoverSV["install_dir"].endswith("share/python3/liftoverSV/bin"):
        g_liftoverSV["install_dir"] = g_liftoverSV["install_dir"][:-29]
    # $LIFTOVERSV/bin/liftoverSV.py
    if g_liftoverSV["install_dir"].endswith("bin"):
        g_liftoverSV["install_dir"] = g_liftoverSV["install_dir"][:-4]

    g_liftoverSV["etc_dir"]  = os.path.join(g_liftoverSV["install_dir"], "etc", "liftoverSV")
    g_liftoverSV["doc_dir"]  = os.path.join(g_liftoverSV["install_dir"], "share", "doc", "liftoverSV")

    pattern = os.path.join(g_liftoverSV["install_dir"], "share", "python*", "liftoverSV")
    matches = glob.glob(pattern)
    python_dir = None
    for dir in matches:
        config_file = os.path.join(dir, "workflow/config.py")
        if os.path.isfile(config_file):
            python_dir = dir
            break
    if python_dir is None:
        raise FileNotFoundError("Cannot find g_liftoverSV['python_dir']")
    g_liftoverSV['python_dir'] = python_dir

    g_liftoverSV["bash_dir"] = os.path.join(g_liftoverSV["install_dir"], "share", "bash", "liftoverSV")

    # Add the correct relative path to sys.path so that Python can locate liftoverSV_*.py
    # The path is computed relative to the location of this script (__file__), ensuring
    # it works regardless of the current working directory from which the script is executed.
    #########################################################################################
    sys.path.insert(0, os.path.join(g_liftoverSV['python_dir']))

    # Import the different modules
    # (to keep here after the definition of the correct relative path to sys.path)
    ##############################################################################
    from workflow.config import configure_liftover_sv
    from core.dna_checks import is_multi_allelic, check_ref_fasta_seq, retrieve_chrom_size
    from workflow.liftover_process import write_the_lifted_vcf

    # Search for the liftoverSV VERSION
    ###################################
    if "version" not in g_liftoverSV:
        runFile = os.path.join(g_liftoverSV["install_dir"], "bin", "liftoverSV.py")
        if os.path.exists(runFile):
            try:
                # Open the file in text mode
                with open(runFile, "rt") as f:
                    for line in f:
                        m = re.match(r"^liftoverSV ([0-9]+\.[0-9]+(\.[0-9]+)?(_beta)?)", line)
                        if m:
                            g_liftoverSV["version"] = m.group(1)
                            break
            except Exception as e:
                print(f"[WARNING] Could not read {runFile}: {e}")


    if "version" not in g_liftoverSV:
        g_liftoverSV["version"] = "X.X"


    # Downloading configuration
    ###########################
    configure_liftover_sv(argv, g_liftoverSV)

    # Display
	#########
    print(f"\nliftoverSV {g_liftoverSV['version']}")
    print("Copyright (C) 2024-current GEOFFROY Veronique")
    print("Please feel free to create a Github issue for any suggestions or bug reports")
    print("https://github.com/lgmgeo/liftoverSV/issues\n")
    print("\nPython version:", platform.python_version(), "\n")
    print("Application name used:")
    print(g_liftoverSV["install_dir"], "\n")



    # Arguments display
    ###################
    print(f"\n[{time.strftime('%H:%M:%S')}] Listing arguments")
    print("           *********************************************")
    print("           liftoverSV has been run with these arguments:")
    print("           *********************************************")

    for key in sorted(g_liftoverSV.keys()):
        if key in ["bash_dir", "doc_dir", "etc_dir", "install_dir", "python_dir", "tcl_dir", "version"]:
            continue
        val = g_liftoverSV[key]
        if val == "":
            continue
        key = key.replace("_", "-")
        print(f"           --{key} {val}")

    print("           *********************************************")


    # Check if the input VCF file contains multi-allelic lines
	##########################################################
    is_multi_allelic(g_liftoverSV)

	# Check if the FASTA REF FILE seems correct (with or without "chr")
	###################################################################
    check_ref_fasta_seq(g_liftoverSV)

	# Memorize the size of the chromosomes in the target build
	#########################################################
    retrieve_chrom_size(g_liftoverSV)    

	# Write the lifted VCF
	######################
    write_the_lifted_vcf(g_liftoverSV)

	# Finished
	##########
    print(f"[{time.strftime('%H:%M:%S')}] Liftover completed successfully.")


if __name__ == "__main__":
    main(sys.argv[1:])



