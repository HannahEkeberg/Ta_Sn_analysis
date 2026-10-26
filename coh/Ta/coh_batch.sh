#!/bin/bash

for D in ./*; do
    if [ -d "$D" ]; then
        cd "$D"
        ./run_coh.sh 
        cd ..
    fi
done
