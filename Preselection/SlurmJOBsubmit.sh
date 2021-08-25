#!/bin/bash
sample=$1
polarity=$2
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G

source  /home/hep/buonaura/conda_setup.sh
conda activate lhcb_ana

#python /home/hep/buonaura/Analysis/RLc/scripts/Preselection/CreateTreePreselectionMC.py --MCTrackerOnly $1 $2 --HLT1_1 --HLT1_2 --presel
python /home/hep/buonaura/Analysis/RLc/scripts/Preselection/CreateTreePreselectionMC.py --MCTrackerOnly $1 $2 --presel

conda deactivate

