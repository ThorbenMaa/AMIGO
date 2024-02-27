#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cp run_AMIGO.sh "$i" # copy water.txt into it
        cp AMIGO_non_interactive.py "$i" # copy water.txt into it
    fi
done