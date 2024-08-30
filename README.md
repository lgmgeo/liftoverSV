<p align="center">
<div align="center">
    <h1 style="font-weight: bold">liftoverSV:
      <h3>Lifts over a Structural Variation VCF file from one reference build to another</h3>
    </h1>
</div>

<br />


## COMMAND LINE USAGE

       $LIFTOVERSV/bin/liftoverSV -I $INPUT_FILE -C $CHAIN_FILE -O $OUTPUT_FILE


## OPTIONS
```bash
--BCFTOOLS,-F <File>          The bcftools path
                              See https://samtools.github.io/bcftools/howtos/install.html
                              Default: "bcftools"

--BEDTOOLS,-B <File>          The bedtools path
                              See https://bedtools.readthedocs.io/en/latest/content/installation.html
                              Default: "bedtools"

--CHAIN,-C <File>             The liftover chain file
                              See https://genome.ucsc.edu/goldenPath/help/chain.html for a description of chain files
                              See http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files.
                              Required

--help,-h <Boolean>           Display the help message
                              Default value: false. Possible values: {true, false}

--INPUTFILE,-I <File>         The SV VCF input file
                              Required

--LIFTOVER,-L <File>          The UCSC Liftover tool path
                              Default: "liftOver"

--OUTPUTFILE,-O <File>        The liftover SV VCF output file
                              Required

--PERCENT,-P <float>          Variation in length authorized for a lifted SV (e.g. difference max between SVLEN < 5%)
                              Default value: 0.05

--REFFASTASEQ,-R <File>       The reference sequence (fasta) for the TARGET genome build (i.e., the new one after the liftover)
```

## How to cite?
Please cite the following doi if you are using this tool in your research:</br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12799803.svg)](https://doi.org/10.5281/zenodo.12799803)

## liftoverSV: 

Lifts over a SV VCF file from one reference build to another:

* Lifts over #CHROM and POS

* Lifts over INFO/END and INFO/SVEND</br>

* Lifts over INFO/SVLEN, INFO/SVSIZE:
   - Lifts over for deletion, duplication and inversion (SVLEN_lifted = End_lifted - Start_Lifted)
   - Keep the same SVLEN/SVSIZE for insertion (the number of the inserted bases remains the same)
   - Set SVLEN/SVSIZE to "." for SV type not equal to DEL, DUP, INV or INS (TRA, CPX...)

* Drop the SV if:</br>
   - Case1: one position (start or end) is lifted while the other doesn't
   - Case2: one position (start or end) goes to a different chrom from the other
   - Case3: "lifted start" > "lifted end"
   - Case4: the distance between the two lifted positions changes significantly (Default: difference between both SVLENs > 5%)

* Update/create and sort some VCF header lines:
	- Checks that the "contig" field includes all the ID attributes (do not include additional optional attributes)</br>
	  e.g. `##contig=<ID=chr22>` added after a lift from chr1 to chr22
	- Create/update the "assembly" field 
	  e.g. `##assembly=liftoverSV used with hg19ToHg38.over.chain`
	- Create the "liftoverSV_version" field
	  e.g. `##liftoverSV_version=0.1.2_beta; hg19ToHg38.over.chain; August 30 2024 12:30`
	- Update the "INFO" and "FORMAT field if one value is missing.</br>
	  As the format (Number, String) is not known, "Number=." and "Type=String" values are used by default:</br>
	  e.g. `##FORMAT=<ID=XXX,Number=.,Type=String,Description="XXX">`</br>
	  e.g. `##INFO=<ID=YYY,Number=.,Type=String,Description="YYY">`

## Requirements
### a) The UCSC Liftover tool (required)
The UCSC Liftover tool needs to be locally installed.</br>
https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
### b) bcftools
The “bcftools” toolset needs to be locally installed to sort the VCF output file
### c) bedtools (to be required in future development)
The “bedtools” toolset will need to be locally installed to lift over sequences (e.g. ACGGTTG]chr1:12569863])

## SV VCF format: Documentation
cf https://samtools.github.io/hts-specs/VCFv4.4.pdf</br>
See section 3: "INFO keys used for structural variants"

## Feature requested for future release

* Checks INFO/CIPOS and INFO/CIEND</br>
=> POS-CIPOS >0</br>
=> END+CIEND < chrom_length</br>

* Lifts over INFO/MEINFO and INFO/METRANS

* Lifts over INFO/HOMLEN, INFO/HOMSEQ
  
* Lifts over ALT when described with square bracket notation. For example, G]17:198982] or ]chr1:3000]A</br>
cf [https://github.com/EUCANCan/variant-extractor](variantextractor) for notation rules



