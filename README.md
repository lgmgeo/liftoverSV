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

## liftoverSV: version 0.1.1_beta

* Lift over #CHROM and POS

* Lift over INFO/END and INFO/SVEND</br>
=> drop the SV if:</br>
   - Case1: one position (start or end) is lifted while the other doesn't
   - Case2: one position (start or end) goes to a different chrom from the other
   - Case3: "lifted start" > "lifted end"
   - Case4: the distance between the two lifted positions changes significantly (difference between both SVLENs > 5%)

* Lift over INFO/SVLEN, INFO/SVSIZE (for deletion, duplication, insertion and inversion)

* The structured contig field includes all the ID attributes (do not include additional optional attributes)</br>
e.g. `##contig=<ID=chr22>`

## Requirements
### a) The UCSC Liftover tool (required)
The UCSC Liftover tool needs to be locally installed.</br>
https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
### b) bedtools (to be required in future development)
The “bedtools” toolset (developed by Quinlan AR) will need to be locally installed to lift over sequences (e.g. ACGGTTG]chr1:12569863])

## SV VCF format: Documentation
cf https://samtools.github.io/hts-specs/VCFv4.4.pdf</br>
See section 3: "INFO keys used for structural variants"

## Feature requested for future release

* Check INFO/CIPOS and INFO/CIEND</br>
=> POS-CIPOS >0</br>
=> END+CIEND < chrom_length</br>

* Liftover INFO/MEINFO and INFO/METRANS

* Liftover INFO/HOMLEN, INFO/HOMSEQ
  
* Liftover ALT when described with square bracket notation. For example, G]17:198982] or ]chr1:3000]A</br>
cf variantextractor for notation rules

