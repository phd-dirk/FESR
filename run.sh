#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`
OUTPUT_FILE="$1/fits.csv"

echo $OUTPUT_FILE
echo -e "$DATE \n"'status,alpha,alphaErr,aGGInv,aGGInvErr,c6,c6Err,c8,c8Err,chi,dof,chiDof,edm' > $OUTPUT_FILE
for i in $( ls $1 ); do
    if [[ $i == *.json ]]; then
        echo $i
        ./build/FESR "$1/$i"
    fi
done
