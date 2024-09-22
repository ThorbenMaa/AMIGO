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

echo start non-PCS data sets
# non-PCS data sets
## declare an array variable
#declare -a arr=("GTB_without_PCS" "MNV_Pdomain_without_PCS" "atcase_fix" "EIN" "HuNoV_Pdomain" "LmUGP_apo_without_PCS" "mbp" "mrsB" "msg" "ubiquitin" "a7a7")
declare -a arr=("ubiquitin")

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
#declare -a arr=("GTB" "LmUGP_with_PCS" "MNV_Pdomain_with_PCS")
declare -a arr=("MNV_Pdomain_with_PCS")

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

cd ..

# synthetic data sets
echo start synthetic MILVproS
cd synthetic_NOE_networks_MILVproS_7_A_2_1ratio


## declare an array variable
#declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")
declare -a arr=("hu_glu_cyclase_apo_chainA")

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
cd ..



echo start synthetic MILVAproS
cd synthetic_NOE_networks_MILVAproS_7A_2_1ratio


## declare an array variable
#declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")
declare -a arr=("hu_glu_cyclase_apo_chainA")

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
cd ..

echo start synthetic MILVATproS
cd synthetic_NOE_networks_MILVATproS_7A_2_1ratio


## declare an array variable
#declare -a arr=("hu_glu_cyclase_apo_chainA" "hu_TYRP1" "hu_TEAD1_chainA" "hu_carbonic_anhydrase_XIII" "hu_ELOVL_fatty_acid_elongase_chainA" "tetanus_toxin_c" "trypsin" "phosphorylase")
declare -a arr=("hu_glu_cyclase_apo_chainA")

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
cd ..


echo fertig