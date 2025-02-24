#!/bin/bash 

set -eo pipefail


CHAIN=$1
REFFASTASEQ=$2


# Check a VCF input file with ALT described with "square bracketed notation"
############################################################################
#
# input not sorted !
# => Ok, output sorted
#
# Manual check of the input/output files:
# => with the UCSC genome browser (e.g. hg19 and "chr13:52465906-52465906" (1-based coordinates)  => C)

rm -f ./output/output_hg38.*

$LIFTOVERSV/bin/liftoverSV -I ./input/input_hg19.vcf -O ./output/output_hg38.vcf -C $CHAIN -R $REFFASTASEQ


gunzip ./output/output_hg38.sort.vcf.gz
if [ -e ./validated_output/validated_output_hg38.sort.vcf.gz ]
then
        gunzip ./validated_output/validated_output_hg38.sort.vcf.gz
fi

compare=`diff ./output/output_hg38.sort.vcf validated_output/validated_output_hg38.sort.vcf || true`

gzip ./output/output_hg38.sort.vcf
gzip ./validated_output/validated_output_hg38.sort.vcf


if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
else
        echo "ok - Finished"
fi


