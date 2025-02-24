#!/bin/bash 

set -eo pipefail


CHAIN=$1
REFFASTASEQ=$2


# Check a VCF input file with ALT described with an "angle bracketed notation"
##############################################################################
#
# ALT = <DEL:SVSIZE=52:AGGREGATED>
#
# INPUT coordinates:  chr22:16848506-16848558
# OUTPUT coordinates: chr22:16367844-16367896
#
# Checks:
# ../scripts/extractDNAseq.tcl chr22 16848506 16848558 $hg19_REFFASTASEQ   => gaatggaatcatcaacgaatggaatcgaatggaatcatcgtctaatggaatca
# ../scripts/extractDNAseq.tcl chr22 16367844 16367896 $hg38_REFFASTASEQ   => gaatggaatcatcaacgaatggaatcgaatggaatcatcgtctaatggaatca



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


