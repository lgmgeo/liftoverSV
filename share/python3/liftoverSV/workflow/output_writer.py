"""
liftoverSV 0.3.1_beta
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
import gzip
import time
from io_tools.file_utils import open_any_text_file, natural_sort_key, print_flush as print
from io_tools.vcf_sorter import VcfSorter

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
def add_new_header_lines(S_new_INFO, S_new_FORMAT, S_new_FILTER, S_lifted_contigs, input_file, g_liftoverSV):
    """
    Add the updated "VCF input_file header" to gzipped output_file:
    - Adding new INFO, FORMAT, and FILTER lines
    - Updating the contig header list (and adding contig lengths from the target build if needed)
    - Adding information in the L_assembly field about the use of liftoverSV
    - Adding the command line used
    """

    print(f"[{time.strftime('%H:%M:%S')}] Writing header in the {g_liftoverSV['output_file']}")
    print(f"           => Updating (if needed) the INFO, FORMAT and FILTER header lines")
    print(f"           => Updating (if needed) the contigs header lines")
    print(f"           => Adding information about the use of liftoverSV")
    print(f"           => Adding/Updating the reference header line")

    # Memorize header info from the "vcf_header_lines.txt" file
    ###########################################################
    header_file = os.path.join(g_liftoverSV['install_dir'], "share", "doc", "liftoverSV", "vcf_header_lines.txt")
    # Format in the "vcf_header_lines.txt" file:
    # - FILTER line: Key Description
    # - INFO line: Key Number Type Description
    # - FORMAT line: Field Number Type Description
    INFO = {}
    FORMAT = {}
    FILTER = {}
    feature = None
    with open_any_text_file(header_file) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue
                
            if line.startswith("#"):
                if line.startswith("# FILTER"):
                    feature = "FILTER"
                elif line.startswith("# INFO"):
                    feature = "INFO"
                elif line.startswith("# FORMAT"):
                    feature = "FORMAT"
                continue

            split_line = line.split()
            if feature == "FILTER":
                key = split_line[0]
                desc = " ".join(split_line[1:])
                FILTER[key] = {"Description": desc}
            elif feature == "INFO":
                key, number, type_ = split_line[:3] ;# The underscore _ is used to avoid using the Python keyword type as a variable name.
                desc = " ".join(split_line[3:])
                INFO[key] = {"Number": number, "Type": type_, "Description": desc}
            elif feature == "FORMAT":
                key, number, type_ = split_line[:3]
                desc = " ".join(split_line[3:])
                FORMAT[key] = {"Number": number, "Type": type_, "Description": desc}


    # Initialisation of the "L_*" variables (1 "L_" variable for each header line category)
    #######################################################################################
    L_file_format, L_ALT, L_META, L_SAMPLE, L_PEDIGREE, L_others = [], [], [], [], [], []
    L_CHROM = None
    # The one to be completed:
    L_contig, L_assembly, L_reference, L_command, L_FILTER, L_INFO, L_FORMAT = [], [], [], [], [], [], []


    # Create the complete contigs header
    ####################################
    sorted_lifted_contigs = sorted(S_lifted_contigs, key=natural_sort_key)
    for contig in sorted_lifted_contigs:
        # Try to get the contig size from g_liftoverSV["size_chrom_target"]
        size_chrom_target = g_liftoverSV["size_chrom_target"].get(contig)
        # Corrects/Add the length of the contig (from source to target) in the line
        if size_chrom_target:
            line = f"##contig=<ID={contig},length={size_chrom_target}>"
        else:
            line = f"##contig=<ID={contig}>"
        # Add the line
        L_contig.append(line)


    # While reading the input VCF, complete the "L_*" header lines 
    ##############################################################
    with open_any_text_file(input_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                if line.startswith("##fileformat="):
                    L_file_format.append(line)
                elif line.startswith("##assembly="):
                    suffix = " (before using the liftoverSV)"
                    if not line.endswith(suffix):
                        line = line + suffix
                    L_assembly.append(line)
                elif line.startswith("##reference="):
                    ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
                    continue
                elif line.startswith("##contig=<ID"):
                    continue
                elif line.startswith("##FILTER=<ID"):
                    L_FILTER.append(line)
                elif line.startswith("##INFO=<ID"):
                    L_INFO.append(line)
                elif line.startswith("##FORMAT=<ID"):
                    L_FORMAT.append(line)
                elif line.startswith("##ALT=<ID"):
                    L_ALT.append(line)
                elif line.startswith("##META=<ID"):
                    L_META.append(line)
                elif line.startswith("##SAMPLE=<ID"):
                    L_SAMPLE.append(line)
                elif line.startswith("##PEDIGREE=<ID"):
                    L_PEDIGREE.append(line)
                elif line.startswith("#CHROM"):
                    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
                    L_CHROM = line
                else:
                    L_others.append(line)
            else:
                break


    # Complete L_FILTER, L_INFO and L_FORMAT (with the given arguments)   
    ###################################################################
    for key in S_new_FILTER:
        if key in FILTER:
            L_FILTER.append(f"##FILTER=<ID={key},Description=\"{FILTER[key]['Description']}\">")
        else:
            L_FILTER.append(f"##FILTER=<ID={key},Description=\"{key}\">")

    for key in S_new_INFO:
        if key in INFO:
            L_INFO.append(f"##INFO=<ID={key},Number={INFO[key]['Number']},Type={INFO[key]['Type']},Description=\"{INFO[key]['Description']}\">")
        else:
            L_INFO.append(f"##INFO=<ID={key},Number=.,Type=String,Description=\"{key}\">")

    for key in S_new_FORMAT:
        if key in FORMAT:
            L_FORMAT.append(f"##FORMAT=<ID={key},Number={FORMAT[key]['Number']},Type={FORMAT[key]['Type']},Description=\"{FORMAT[key]['Description']}\">")
        else:
            L_FORMAT.append(f"##FORMAT=<ID={key},Number=.,Type=String,Description=\"{key}\">")

    # Create the "L_reference" and "L_command"
    ##########################################
    L_reference = [f"##reference={g_liftoverSV['ref_fasta_seq']}"]
    
    ##command
    L_command.append(f"##liftoverSV_command={sys.executable} " + " ".join(sys.argv))
    L_command.append(f"##liftoverSV_version={g_liftoverSV['version']}")


    # Combine all header lines
    ##########################
    L_to_write = []
    for group in [L_file_format, L_reference, L_assembly, L_command, L_contig, L_FILTER, L_INFO, L_FORMAT, L_ALT, L_META, L_SAMPLE, L_PEDIGREE, L_others]:
        if group:
            L_to_write.extend(sorted(group, key=natural_sort_key))
    
    if L_CHROM is None:
        print("\nERROR: Header line starting with #CHROM not found in the VCF!")
        print("   This file is not a valid VCF according to the specification.")
        print("   liftoverSV cannot continue.\n")
        sys.exit(1)
    else:
        L_to_write.append(L_CHROM)

    # Write header to gzipped output_file
    #####################################
    updated_header="\n".join(L_to_write)  
    del L_to_write

    # Remove the output_file if already exists
    if os.path.exists(g_liftoverSV['output_file']):
        try:
            os.remove(g_liftoverSV['output_file'])
        except OSError as e:
            print(f"Warning: cannot remove {g_liftoverSV['output_file']}: {e}")

    # Write header to gzipped output_file
    with gzip.open(g_liftoverSV['output_file'], "wt") as out_f:
        # write header 
        out_f.write(updated_header + "\n")



def sort_and_compress_the_lifted_vcf(tmp_output_file, g_liftoverSV):
    """
    Do not modify the header in g_liftoverSV['output_file'].
    Sort and compress the tmp output VCF file and write variant lines to g_liftoverSV['output_file']
    """

    # Sort and compress the output file 
    sorter = VcfSorter(
        vcf_to_sort=tmp_output_file,
        sorted_vcf=g_liftoverSV['output_file'],
        overwrite=False
    )
    sorter.sort(g_liftoverSV)

    # Clean: remove tmp output if output_file exists
    if os.path.exists(g_liftoverSV['output_file']):
        try:
            print(f"[{time.strftime('%H:%M:%S')}] Removing {tmp_output_file}")
            os.remove(tmp_output_file)
        except OSError as e:
            print(f"Warning: cannot remove {tmp_output_file}: {e}")