#! /bin/bash
### Submit this Script with: sbatch script.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=longterm
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 16
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=14-00:00:00
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=10GB
#  Find your job easier with a name:
#SBATCH --job-name=GTB

echo start non-PCS data sets
# non-PCS data sets
## declare an array variable
declare -a arr=("atcase_fix" "EIN" "HuNoV_Pdomain" "LmUGP_apo_without_PCS" "mbp" "mrsB" "msg" "ubiquitin" "a7a7" )

## now loop through the above array
for i in "${arr[@]}"
do
    if [ -d "$i" ] # if it's a directory
    then
        echo testing $i now...
        
        cd $i
        cd popt
        bash sbatch_test.sh
        cd ..
        cd p\=1
        bash sbatch_test.sh
        cd ..
        cd ..
        
    fi
done

# PCS data sets

echo start PCSs data sets
cd PCS_data_sets
## declare an array variable
declare -a arr=("GTB" "LmUGP_with_PCS" "MNV_Pdomain_with_PCS")

## now loop through the above array
for i in "${arr[@]}"
do
    if [ -d "$i" ] # if it's a directory
    then
        echo testing $i now...
        
        cd $i
        bash sbatch_test.sh
        cd ..
        
    fi
done


echo fertig