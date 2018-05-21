#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`

mkdir ./output/$DATE
 echo -e 'numS,alpha,alphaErr,aGGInv,aGGInvErr,rhoVpA,rhoVpAErr,cVpA,cVpAErr,chi,chiDof,edm,sSet,weight' > ./output/$DATE/fits.csv
 for i in $( ls ./configurations/ ); do
     # run fit
     cp ./configurations/$i ./configuration.json
     ./build/FESR ./output/$DATE/fits.csv
     # output
     cp ./configurations/$i ./output/$DATE/$i
 done
