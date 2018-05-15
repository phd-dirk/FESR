#!/bin/bash
DATE=`date '+%Y-%m-%d_%H:%M:%S'`

mkdir ./output/$DATE
 echo -e 'alpha \t alpha_error \t aGGInv \t aGGinv_error \t rhoVpA \t rhoVpA_error \t c8VpA \t c8VpA_error \t chi2 \t edm \t s0s \t weight' > ./output/$DATE/fits.dat
 for i in $( ls ./configurations/ ); do
     # run fit
     cp ./configurations/$i ./configuration.json
     ./build/FESR ./output/$DATE/fits.dat
     # output
     cp ./configurations/$i ./output/$DATE/$i
 done
