#!/bin/bash

# clean
rm ./output/*
echo -e 'alpha \t alpha_error \t aGGInv \t aGGinv_error \t rhoVpA \t rhoVpA_error \t c8VpA \t c8VpA_error' > ./output/fits.dat
for i in $( ls ./configurations/ ); do
    # run fit
    cp ./configurations/$i ./configuration.json
    ./build/FESR
    # output
    cp ./configurations/$i ./output/$i
done
