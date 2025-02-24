#!/bin/bash 

set -eo pipefail


CHAIN=$1
REFFASTASEQ=$2


rm -f ./output/output_hg38.*

# Check a compressed VCF in input
$LIFTOVERSV/bin/liftoverSV -I ./input/input_hg19.vcf.gz -O ./output/output_hg38.vcf.gz -C $CHAIN -R $REFFASTASEQ

# Check different equivalent output file names
$LIFTOVERSV/bin/liftoverSV -I ./input/input_hg19.vcf.gz -O ./output/output_hg38.vcf.gz -C $CHAIN -R $REFFASTASEQ
$LIFTOVERSV/bin/liftoverSV -I ./input/input_hg19.vcf.gz -O ./output/output_hg38.sort.vcf.gz -C $CHAIN -R $REFFASTASEQ



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

