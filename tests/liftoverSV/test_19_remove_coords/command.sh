#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2

# INPUT:
########

# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
# chr17   198982  trn_no_mateid_a G       G]chr2:321681]  .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1
# chr2    321681  trn_no_mateid_b G       G]chr17:198982] .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1
# chr13   53040041        ins_by_gridss   T       TATATATATACACAC[chr13:53040042[ .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1
# chr13   53040042        ins_by_gridss   A       ]chr13:53040041]ATATATATACACACA .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1
# chr15   53040041        empty_del_a     A       A[chr15:53040042[       .       PASS    SVTYPE=BND;EXTRA=DEL_NO_SIZE;MATEID=empty_del_b GT      0/0     0/1
# chr15   53040042        empty_del_b     T       ]chr15:53040041]T       .       PASS    SVTYPE=BND;EXTRA=DEL_NO_SIZE;MATEID=empty_del_a GT      0/0     0/1


rm -f ./output/output_hg38.*

python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.1.vcf -o ./output/output_hg38.1.vcf -c $chain -r $ref_fasta_seq -R

gunzip ./output/output_hg38.1.sort.vcf.gz
if [ -e ./validated_output/validated_output_hg38.1.sort.vcf.gz ]
then
        gunzip ./validated_output/validated_output_hg38.1.sort.vcf.gz
fi

compare=`diff ./output/output_hg38.1.sort.vcf validated_output/validated_output_hg38.1.sort.vcf || true`

gzip ./output/output_hg38.1.sort.vcf
gzip ./validated_output/validated_output_hg38.1.sort.vcf


if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values (input_hg19.1.vcf)"
fi


python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.2.vcf -o ./output/output_hg38.2.vcf -c $chain -r $ref_fasta_seq -R -D "SVLEN"

gunzip ./output/output_hg38.2.sort.vcf.gz
if [ -e ./validated_output/validated_output_hg38.2.sort.vcf.gz ]
then
        gunzip ./validated_output/validated_output_hg38.2.sort.vcf.gz
fi

compare=`diff ./output/output_hg38.2.sort.vcf validated_output/validated_output_hg38.2.sort.vcf || true`

gzip ./output/output_hg38.2.sort.vcf
gzip ./validated_output/validated_output_hg38.2.sort.vcf


if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values (input_hg19.2.vcf)"
fi
echo "ok - Finished"
