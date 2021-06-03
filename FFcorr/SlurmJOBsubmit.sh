#!/bin/bash
sample=$1
polarity=$2
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

#source  /home/hep/buonaura/conda_setup.sh
#conda activate lhcb_ana
source ~/myroot/bin/thisroot.sh

python3 /home/hep/buonaura/Analysis/RLc/scripts/FFcorr/Store_FFcorr.py --MCTrackerOnly $1 $2 

#conda deactivate

