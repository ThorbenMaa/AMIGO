#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cp AMIGO_non_interactive.py "$i" # copy water.txt into it
        cp 2oef.pdb "$i"

    fi
done