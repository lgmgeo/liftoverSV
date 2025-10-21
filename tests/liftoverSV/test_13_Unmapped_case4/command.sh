#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2


rm -f ./output/output_hg38.*

python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf -o ./output/output_hg38.vcf -c $chain -r $ref_fasta_seq

compare=`diff ./output/output_hg38.unmapped validated_output/validated_output_hg38.unmapped || true`


if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
else
        echo "ok - Finished"
fi

