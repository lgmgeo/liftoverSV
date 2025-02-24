#!/bin/bash 

set -eo pipefail


CHAIN=$1
REFFASTASEQ=$2

# INPUT: => 1 SV (1 unmapped)
# chr10   46951150        .       TGTGTGTGTGTGTGTGTGTGTG  [chrX:96440306[ATATAT   70950   GraphtyperFilter;inf50bp        END=46951170;AN=1698;AC=900     GT      0/1     0/1

# unmapped INPUT:
# CHROM=chr10
# POS=46951150
# END=46951170
# => POS < END
#
# Lift of the input coordinates:
# CHROM=chr10
# POS=46598467
# END=46598447
# => POS > END

# => Pas de SV dans l'output !!!
# => Affiche un message dans la sortie de liftoverSV

rm -f ./output/output_hg38.*

$LIFTOVERSV/bin/liftoverSV -I ./input/input_hg19.vcf -O ./output/output_hg38.vcf -C $CHAIN -R $REFFASTASEQ > liftoverSV.log

result=`grep -c "...No SV lifted!" liftoverSV.log || true`

 
if [ "$result" ]
then
        echo "ok - Finished"
else
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
fi


