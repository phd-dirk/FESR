#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`
OUTPUT_FILE="$1/fits.csv"

echo $OUTPUT_FILE
echo "# $DATE" > $OUTPUT_FILE
echo "status,alpha,alphaErr,aGGInv,aGGInvErr,c6,c6Err,c8,c8Err,c10,c10Err,c12,c12Err,chi,dof,chiDof,edm,del^(0),del^(4),del^(6),del^(8),del^(10),del^(12)" >> $OUTPUT_FILE
for i in $( ls $1 ); do
    if [[ $i == *.json ]]; then
        echo $i
        nice ./build/FESR "$1/$i"
    fi
done
