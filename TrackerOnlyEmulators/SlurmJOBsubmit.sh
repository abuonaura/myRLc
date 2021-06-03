#!/bin/bash
#SBATCH --job-name=CreateTriggerMCTO_files_%spl_%pol
#SBATCH --output=/home/hep/buonaura/Analysis/RLc/scripts/TrackerOnlyEmulators/logs/log_%spl_%pol.txt
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
source  /home/hep/buonaura/conda_setup.sh
conda activate lhcb_ana

python /home/hep/buonaura/Analysis/RLc/scripts/TrackerOnlyEmulators/run_triggers_emulators.py -f $1 --L0HTOS --L0GTIS --HLT1

conda deactivate

