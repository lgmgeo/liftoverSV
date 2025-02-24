#!/bin/bash

CHAIN=$1
REFFASTASEQ=$2

for f in test_*/command.sh
do
    if [[ "$f" == *00* ]]
    then
        continue
    fi
    d=`dirname $f`
    cd $d
    echo "######## $d "
    echo "###############################################################################"
    date
    ./command.sh $CHAIN $REFFASTASEQ >& command.log
    grep -i error command.log | grep -v "Exit without error"

    # diff
    egrep "(^> ##)|(^< ##)|(---)" command.log

    grep "ok - Finished" command.log

    echo ""
    cd ..
done


echo "ok - Finished all checkings"

