#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2

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

python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf -o ./output/output_hg38.vcf -c $chain -r $ref_fasta_seq > liftoverSV.log

result=`grep -c "* 0 mapped SV" liftoverSV.log || true`

 
if [ "$result" ]
then
        echo "ok - Finished"
else
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
fi


