#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        cd $i
        sbatch sbatch_test.sh
        cd ..
    fi
done
echo fertig