#!/bin/bash
for i in ./* # iterate over all files in current dir
do
    if [ -d "$i" ] # if it's a directory
    then
        echo $i
        
        cd $i
        cd popt
        sbatch sbatch_test.sh
        cd ..
        cd p\=1
        sbatch sbatch_test.sh
        cd ..
        cd ..
        
    fi
done
echo fertig