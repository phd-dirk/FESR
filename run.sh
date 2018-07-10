#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`
OUTPUT_DIR="$DATE""_$1"
echo $OUTPUT_DIR

mkdir ./output/$OUTPUT_DIR
echo -e 'numS,alpha,alphaErr,aGGInv,aGGInvErr,rhoVpA,rhoVpAErr,cVpA,cVpAErr,chi,chiDof,edm,sSet,weight' > ./output/$OUTPUT_DIR/fits.csv
for i in $( ls ./configurations/ ); do
    # run fit
    cp ./configurations/$i ./configuration.json
    ./build/FESR ./output/$OUTPUT_DIR/fits.csvdh
    # output
    cp ./configurations/$i ./output/$OUTPUT_DIR/$i
done
