#! /bin/bash

# popt folders
declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")

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


# p1 folders
declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")

## now loop through the above array
for i in "${arr[@]}"
do
    if [ -d "$i" ] # if it's a directory
    then
        echo generating $i now...
        
        for j in 7-5 8 8-5
        do
          #cd $j
          bash simulate.sh $i/p\=1/$j *.pdb
        done
    fi
done