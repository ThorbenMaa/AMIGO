#! /bin/bash

# popt folders
declare -a arr=("phosphorylase")

## now loop through the above array
for i in "${arr[@]}"
do
    if [ -d "$i" ] # if it's a directory
    then
        echo generating $i now...
        
        for j in 7 7-5 8 8-5 9 9-5 10
        do
          #cd $j
          bash simulate.sh $i/popt/$j *.pdb
        done
    fi
done


