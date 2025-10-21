#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2


rm -f ./output/output_hg38.*

# Check a compressed VCF in input

# Check different equivalent output file names
python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf.gz -o ./output/output_hg38.vcf.gz -c $chain -r $ref_fasta_seq
python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf.gz -o ./output/output_hg38.sort.vcf.gz -c $chain -r $ref_fasta_seq
python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf.gz -o ./output/output_hg38 -c $chain -r $ref_fasta_seq



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

