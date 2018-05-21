#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`

mkdir ./output/$DATE
 echo -e '#s0,alpha,alpha_error,aGGInv,aGGinv_error,rhoVpA,rhoVpA_error,c8VpA,c8VpA_error,chi2,chi2/dos,edm,s0s,weight' > ./output/$DATE/fits.csv
 for i in $( ls ./configurations/ ); do
     # run fit
     cp ./configurations/$i ./configuration.json
     ./build/FESR ./output/$DATE/fits.csv
     # output
     cp ./configurations/$i ./output/$DATE/$i
 done
