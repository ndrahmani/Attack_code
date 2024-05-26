#!/bin/bash

#SBATCH --job-name=Attack_demonstration

#SBATCH -o %x-%j.out

#SBATCH -e %x-%j.err

#SBATCH --partition=mediumq
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # 1 task per node
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G               # Memory per CPU
#SBATCH --time=0-15:00:00       # Estimated completion time





module load Anaconda3
export CONDA_ENVS_PATH=/data/$USER/envs
source activate /data/n.rahmani/envs/sageg6k

gp -q attack_demonstration.gp