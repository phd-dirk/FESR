#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`
OUTPUT_FILE="$1/fits.csv"

echo $OUTPUT_FILE
echo -e "$DATE \n"'status,alpha,alphaErr,aGGInv,aGGInvErr,O6,O6Err,O8,O8Err,chi,dof,chiDof,edm' > $OUTPUT_FILE
for i in $( ls $1 ); do
    if [ "$i" != "fits.csv" ]; then
        echo $i
        ./build/FESR "$1/$i" "$1/fits.csv"
    fi
done
