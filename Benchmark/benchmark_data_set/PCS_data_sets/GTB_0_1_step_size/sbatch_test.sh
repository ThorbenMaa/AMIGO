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
#SBATCH --time=7-00:00:00
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=10GB
#  Find your job easier with a name:
#SBATCH --job-name=GTB


# Initialize the module system:
source /etc/profile.d/modules.sh


# ADJUST DIRECTORY HERE!!!
echo start

echo running on a
lscpu

start=`date +%s`
echo started $start

cd ./6-5
bash run_AMIGO.sh
cd ./../7
bash run_AMIGO.sh
cd ./../7-5
bash run_AMIGO.sh
#cd ./../9-5
#bash run_AMIGO.sh
#cd ./../10
#bash run_AMIGO.sh


end=`date +%s`
echo ended $end

runtime=$((end-start))
echo runtime $runtime

echo fertig 