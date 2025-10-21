#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2

# Single breakend:
# Definition of a breakend that is not part of a novel adjacency. 
# (called "single breakends", because they lack a mate)
# e.g. REF=G ALT=.G (ou ALT=G.)

# If an insertion is detected but only the first few base-pairs provided by overhanging reads could be
# assembled, then this inserted sequence can be provided on that line.
# e.g. REF=A ALT=.TGCA


rm -f ./output/output_hg38.*

python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf -o ./output/output_hg38.vcf -c $chain -r $ref_fasta_seq


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


