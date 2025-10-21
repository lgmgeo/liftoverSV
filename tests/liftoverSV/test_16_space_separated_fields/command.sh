#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2


python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf -o ./output/output_hg38.vcf -c $chain -r $ref_fasta_seq > liftoverSV.log || true

result=`grep -c "Incorrect number of fields" liftoverSV.log || true`


if [ "$result" ]
then
        echo "ok - Finished"
else
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
fi
  


