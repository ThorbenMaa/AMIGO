#! /bin/bash
### Submit this Script with: sbatch script.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 16
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=2-00:00:00
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=10GB
#  Find your job easier with a name:
#SBATCH --job-name=GTB


# Initialize the module system:
source /etc/profile.d/modules.sh


# ADJUST DIRECTORY HERE!!!
echo start

cd ./8
bash run_AMIGO.sh
cd ./../8-5
bash run_AMIGO.sh
cd ./../9-5
bash run_AMIGO.sh





echo fertig 