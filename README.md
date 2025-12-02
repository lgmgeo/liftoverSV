
<div align="center">
  <h1 style="font-weight: bold; margin-bottom: 0.2em;">liftoverSV:</h1>
  <h3 style="margin-top: 0;">Lifts over a Structural Variation VCF file from one QUERY reference build to a TARGET reference build</h3>
</div>




- [Requirements](#requirements)
- [Quick Installation](#quick-installation)
- [Command line usage / Options](#command-line-usage--options)
- [Outputs](#outputs)
- [How to cite?](#how-to-cite)
- [Coordinate, sequence and metadata management](#coordinate-sequence-and-metadata-management)
- [Criteria for dropping SVs during liftover](#criteria-for-dropping-svs-during-liftover)
- [Handling IDs for lifted SVs](#handling-ids-for-lifted-svs)
- [Examples of authorized formats](#examples-of-authorized-formats)
- [Tests](#tests)
- [SV representation in VCF files](#sv-representation-in-vcf-files)


## Requirements

liftoverSV requires the following Python packages:
```
pyfaidx==0.9.0.3
pyliftover==0.4.1
```
These dependencies ensure efficient VCF parsing, sequence extraction, liftover operations, and high-performance data handling.



## Quick Installation

The sources can be cloned to any directory:
```
cd /path/to/install/
git clone git@github.com:lgmgeo/liftoverSV.git

```


## Command line usage / Options

```bash
usage: liftoverSV.py [-h] [-V] -c <File> -i <File> -r <File> [-d <Dir>] -o <File> [-w N_WORKERS] [-z CHUNK_SIZE] [-p <float>] [-v]


optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit

Input files:
  -c <File>, --chain <File>
                        the liftover chain file
                        see https://genome.ucsc.edu/goldenPath/help/chain.html for a description of chain files
                        see http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files
                        required
  -i <File>, --input-file <File>
                        the SV VCF input file
                        gzipped VCF file is supported
                        multi-allelic lines are not allowed
                        required
  -r <File>, --ref-fasta-seq <File>
                        the reference sequence (fasta) for the TARGET genome build (i.e. the new one after the liftover)
                        required

Output options:
  -d <Dir>, --output-dir <Dir>
                        the liftover SV VCF output directory
                        default: current directory
  -o <File>, --output-base-name <File>
                        Base name for output (generates FILE.sort.vcf.gz and FILE.unmapped)
                        required

Performance:
  -w N_WORKERS, --n-workers N_WORKERS
                        number of parallel worker processes to use for liftover.
                        increasing the number of workers can speed up processing on multi-core machines.
                        default: 8
  -z CHUNK_SIZE, --chunk-size CHUNK_SIZE
                        number of VCF lines to process per chunk.
                        processing the VCF in chunks reduces memory usage and enables parallel liftover.
                        default: 50000

Behavior:
  -p <float>, --percent <float>
                        variation in length authorized for a lifted SV (e.g. difference max between both SVLENs < 5%)
                        default value: 0.05
  -T <Dir>, --tmp-dir <Dir>
                        Directory where temporary files will be created.
                            If not provided, the system default temporary directory is used.
  -v, --verbose         enable verbose output
```

## Outputs
Running the tool will generate two output files:
* File.sort.vcf.gz - the sorted and compressed VCF file containing the successfully lifted SVs
* File.unmapped - a report detailing all SVs that could not be lifted


## How to cite?
Please cite the following doi if you are using this tool in your research:</br>
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.14922213-blue.svg)](https://doi.org/10.5281/zenodo.14922213)<br>


## Coordinate, sequence and metadata management
liftoverSV lifts over a SV VCF file from one QUERY reference build to a TARGET reference build.

* Lift over #CHROM, POS, REF, ALT, INFO/END and INFO/SVEND</br>
   (genomic coordinates and sequences are lifted)

* Lift over INFO/SVLEN, INFO/SVSIZE:
   - Lifts over for deletion, duplication and inversion (SVLEN_lifted = End_lifted - Start_lifted)
   - Keep the same SVLEN/SVSIZE for insertion (the number of the inserted bases remains the same)
   - Set SVLEN/SVSIZE to "." for SV type not equal to DEL, DUP, INV or INS (TRA, CPX...)

* Check/Update INFO/CIPOS and INFO/CIEND, so that with the target reference build:
   - POS-CIPOS >= 1 (VCF coordinates are 1-based)</br>
   - END+CIEND <= chromosome length

* Update/create and sort some VCF header lines:
    - Checks that the "contig" field includes all the ID attributes</br>
      e.g. `##contig=<ID=chr22,length=50818468>` added after a lift from chr1 to chr22
    - Create/update the "reference" field</br>
      e.g. `##reference=hg38.fa
    - Add the liftoverSV command line used</br>
      => useful for reproducibility and traceability
    - Create the "liftoverSV_version" field</br>
      e.g. `##liftoverSV_version=1.0.0
    - Update the "INFO", "FORMAT" and "FILTER" fields if one value is missing.</br>
      Possible rescues in the \"$liftoverSV/share/doc/liftoverSV/vcf_header_lines.txt"\ file.</br>
      Else, as the format (Number, String) is not known, "Number=." and "Type=String" values are used by default:</br>
      e.g. `##FORMAT=<ID=XXX,Number=.,Type=String,Description="XXX">`</br>
      e.g. `##INFO=<ID=YYY,Number=.,Type=String,Description="YYY">`


## Criteria for dropping SVs during liftover
Drop the SV if:

	Case 1: One or more required positions fail to lift:
	- "Start" coordinate
	- "End" coordinate
	- "Last REF base" coordinate
	- "Square-bracketed ALT" coordinate
	- "Last base in REF before a deleted segment" coordinate

	Case 2: Lifted positions map to different chromosomes (except for translocations):
	- "Start" and "end" coordinates
	- "Start" and "last REF base" coordinates
	- "Start" and "square-bracketed ALT" coordinates
	- "Start" and "last base in REF before a deleted segment" coordinates

	Case 3: Reversed order between lifted positions:
	- "Start" and "end" coordinates
	- "Start" and "square-bracketed ALT" coordinates

	Case 4: Significant change in distance after liftover (default: >5% of SVLEN):
	- Distance between "start" and "end" coordinates changed significantly
	- Distance between "start" and "square-bracketed ALT" coordinates changed significantly

	Case 5: Complex or inconsistent REF/ALT sequences
	- REF contains '.' or '*' inside the sequence
	- ALT contains '.' or '*' inside the sequence
	- Deletion: ALT sequence not at the beginning of REF sequence (e.g., REF="ATTCTTG", ALT="TC")
	- Insertion: REF sequence not at the beginning of ALT sequence (e.g., REF="TC", ALT="ATTCTTG")
	- Insertion with single breakend: REF sequence not opposite to '.' in ALT (e.g., REF="G", ALT=".TTTTTTC")
	- Square-bracketed ALT badly formatted
	- The REF sequence differs from the original after liftover (see REF and ALT)                                                                                                                                                                    
All details on dropped SVs are recorded in the `output_file.unmapped` file.


## Handling IDs for lifted SVs
During liftover of a VCF file, some variants may have missing IDs (i.e. ID = ".").

To maintain traceability:<br>
* If a variant already has an ID, the original ID is retained in the output VCF
* If a variant has a missing ID ("."), the pipeline automatically generates a unique ID in the output VCF:<br>
	ID = `lifted_from_<line_number>`<br>
	where `<line_number>` is the line in the input VCF where the original structural variant (SV) appears

This ensures that every lifted variant has a consistent and unique identifier for downstream analyses, while preserving existing IDs when present.

WARNING:<br> 
In some VCF files, the ID field is formatted with a combination of chromosome and position (human-readable ID): <br>
e.g. `chr22:16848506:DG`<br>
Since the ID values are essential to track SV across the different genome build versions, no lift is done on these IDs.


## Examples of authorized formats
| CHROM | POS      | REF                | ALT                             | FILTER | INFO                                              |
| :---: | :------: | :----------------: | :-----------------------------: | :----: | :-----------------------------------------------: |
| chr22 | 16848506 | G                  | &lt;DEL&gt;                     | PASS   | END=16848558;SVLEN=52;SVSIZE=52;SVTYPE=DEL        |
| chr17 | 198982   | G                  | G]chr2:321681]                  | PASS   | SVTYPE=BND;EXTRA=TRA_PAIRED_WITHOUT_MATE_ID       |
| chr2  | 321681   | G                  | G]chr17:198982]                 | PASS   | SVTYPE=BND;EXTRA=TRA_PAIRED_WITHOUT_MATE_ID       |
| chr13 | 53040041 | T                  | TATATATATACACAC[chr13:53040042[ | PASS   | SVTYPE=INS                                        |
| chr1  | 2523792  | ACGCCCCCTCCCCTGCTGTGCTGGCACCC<br/>CCTCCCCTGCCGCGCTGATGCCCCCTCCC<br/>CTGATGCACTGGCGCCCCCTCCCCTGCCA<br/>TGCTGACGCCCCCTCCCCTGCCGTGCTGG<br/>CGCCCCCTCCCC  | A      | PASS | VARTYPE=SV;SVTYPE=DEL;SVLEN=-128 |
| chr1  | 2524045  | A                  | ATGCCCCCTCCCCTGAGGCACTGGTGCCC<br/>CCCTCCCCTGCAGCGCTGATGCCCCCCCTC<br/>CCCTGCCATGCTGACGCCCCCTCCCCTGAT<br>GCACTGG | LowQUAL | VARTYPE=SV;SVTYPE=INS;SVLEN=95 |
| chr1  | 150625563| C                  | .TTTTTTC                        | PASS   | EVENT=single_breakend|INS                         |
| chr2  |  321681  |    G               |  G.                             | PASS   | EVENT=single_breakend


## Tests
A set of tests is included to ensure the correct functioning of liftoverSV. <br>
The test scripts are located in the tests/liftoverSV/ directory. <br> <br>

Running them after installation or before submitting changes is recommended to verify that everything works as expected.


## SV representation in VCF files
cf https://samtools.github.io/hts-specs/VCFv4.4.pdf</br>
See section 3: "INFO keys used for structural variants"




