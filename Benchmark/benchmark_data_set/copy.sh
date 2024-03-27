#!/bin/bash
for i in ./*/p=1* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cp compare.py "$i" # copy water.txt into it
    fi
done