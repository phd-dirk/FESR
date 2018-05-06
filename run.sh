#!/bin/bash
for i in $( ls ./configurations/ ); do
   cp ./configurations/$i ./configuration.json
   ./build/FESR
done
