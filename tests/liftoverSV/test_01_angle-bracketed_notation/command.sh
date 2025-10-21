#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2


# Check a VCF input file with ALT described with an "angle bracketed notation"
##############################################################################
#
# ALT = <DEL:SVSIZE=52:AGGREGATED>
#
# INPUT coordinates:  chr22:16848506-16848558
# OUTPUT coordinates: chr22:16367844-16367896
#
# Checks:
# ../scripts/extractDNAseq.py chr22 16848506 16848558 $hg19_ref_fasta_seq   => gaatggaatcatcaacgaatggaatcgaatggaatcatcgtctaatggaatca
# ../scripts/extractDNAseq.py chr22 16367844 16367896 $hg38_ref_fasta_seq   => gaatggaatcatcaacgaatggaatcgaatggaatcatcgtctaatggaatca

# The input VCF header is badly formatted : 
# Header with "chr1" (##contig=<ID=chr1,length=249250621>) in place of "chr22"



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

echo ""
echo ""


if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
else
        echo "ok - Finished"
fi


