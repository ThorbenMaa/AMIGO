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
# # request special node
# #SBATCH --nodelist=node01
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=10GB
#  Find your job easier with a name:
#SBATCH --job-name=GTB



echo start synthetic MILVATproR



## declare an array variable
declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")


## now loop through the above array
for i in "${arr[@]}"
do
    if [ -d "$i" ] # if it's a directory
    then
        echo testing $i now...
        
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