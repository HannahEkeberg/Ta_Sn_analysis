#!/bin/bash

for D in ./*; do
    if [ -d "$D" ]; then
        cd "$D"
        ls
        runE defaults 
        sed -i 's/FITLEV   1./FITLEV   0./' defaults.inp
        runE defaults 
        cd ..
    fi
done
