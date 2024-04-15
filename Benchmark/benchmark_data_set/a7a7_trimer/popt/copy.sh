#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cp 1yau_single_ring_ABC_chain.pdb "$i" # copy water.txt into it
    fi
done