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


# Size of chr22 with the query build (hg38):
############################################
# 50818468
# (indicated in the "hg19ToHg38.over.chain" file)


# Output VCF
############
# We should obtain in the output:

# chr22   16367844        chr22:16848506:dg       n       <del:svsize=52:aggregated>      1643    pass    end=16367896;svlen=52;svsize=52;svtype=del;cipos=-16367843,0;ciend=0,34450572        gt      0/1     1/1

#
# chrom = chr22
# pos = 16367844
# end = 16367896
# svlen = 52
# svsize = 52
# cipos = -16367843,0 
# ciend = 0,34450572 


CHAIN=$1
REFFASTASEQ=$2


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


