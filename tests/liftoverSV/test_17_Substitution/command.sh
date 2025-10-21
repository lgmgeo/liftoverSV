#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2


# INPUT:
########

# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
# chr5    999888  SNV-seq_ref_ok  C       G       .       PASS    AC=15   GT      0/0     0/1
# chr5    999888  SNV-seq_ref_bad A       G       .       PASS    AC=15   GT      0/0     0/1

# hg19, chr5:999888, "C"
#	=> REF ok for SNV-seq_ref_ok: C
#	=> Bad REF for SNV-seq_ref_bad: A




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


