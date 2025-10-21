#!/bin/bash 

set -eo pipefail


# Input VCF
###########
# chr22   16848506        chr22:16848506:DG       N       <DEL:SVSIZE=52:AGGREGATED>      1643    PASS    END=16848558;SVLEN=52;SVSIZE=52;SVTYPE=DEL;CIPOS=-16500000,0;CIEND=0,34450600        GT      0/1     1/1
#
# chrom = chr22
# POS = 16848506
# END = 16848558
# SVLEN = 52
# SVSIZE = 52
# CIPOS = -16500000,0
# CIEND = 0,34450600


# Size of chr22 with the target build (hg38):
############################################
# 50818468
# (indicated in the "hg19ToHg38.over.chain" file)


# Output VCF
############
# chr22   16367844   chr22:16848506:DG    N    <DEL:SVSIZE=52:AGGREGATED>    1643    PASS    END=16367896;SVLEN=52;SVSIZE=52;SVTYPE=DEL;CIPOS=-16367843,0;CIEND=0,34450572        GT      0/1     1/1
#
# chrom = chr22
# POS = 16367844
# END = 16367896
# SVLEN = 52
# SVSIZE = 52
# CIPOS = -16367843,0 = -(POS - 1)),0
# CIEND = 0,34450572 = 0,(chr22_size - END) = 0,(50818468 - 16367896)

chain=$1
ref_fasta_seq=$2


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


