#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cp noe.txt "$i" # copy water.txt into it
    fi
done