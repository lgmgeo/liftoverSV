#!/bin/bash 

set -eo pipefail


chain=$1
ref_fasta_seq=$2

# INPUT: => 2 SV (1 unmapped + 1 mapped)
# chr10   46951150        .       TGTGTGTGTGTGTGTGTGTGTG  [chrX:96440306[ATATAT   70950   GraphtyperFilter;inf50bp        END=46951170;AN=1698;AC=900     GT      0/1     0/1
# chr22   16848506        .       N       <DEL:SVSIZE=52:AGGREGATED>      1643    PASS    END=16848558;SVLEN=52;SVSIZE=52;SVTYPE=DEL      GT      0/1     1/1

# unmapped SV from the INPUT:
# CHROM=chr10
# POS=46951150
# END=46951170
# => POS < END
#
# Lift of the coordinates:
# CHROM=chr10
# POS=46598467
# END=46598447
# => POS > END


rm -f ./output/output_hg38.*

python3 $LIFTOVERSV/bin/liftoverSV.py -i ./input/input_hg19.vcf -o ./output/output_hg38.vcf -c $chain -r $ref_fasta_seq


compare=`diff ./output/output_hg38.unmapped validated_output/validated_output_hg38.unmapped || true`

if [ "$compare" ]
then
        echo "$compare"
        echo `basename $(pwd)`": ERROR, not the expected values"
else
        echo "ok - Finished"
fi

