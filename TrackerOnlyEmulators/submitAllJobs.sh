#!/bin/bash

samples=("Lb_Lcmunu" "Lb_Lctaunu" "Lb_LcDs" "Lb_Lc2625munu" "Lb_Lc2625taunu" "Lb_Lc2625Ds" "Lb_Lc2593munu" "Lb_Lc2593taunu" "Lb_Lc2593Ds" "B_Lcpbarmunu" "Lb_Lc2765munu" "Lb_Lc2880munu")
polarities=("MagUp" "MagDown")

for sample in "${samples[@]}";
do
	echo "$sample"
	export spl="$sample"
	for polarity in "${polarities[@]}";
	do
		echo "$polarity"
		export pol="$polarity"
		inputFile="/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/"$sample"_"$polarity".root"
		echo "$inputFile"
		#python run_triggers_emulators.py -f "$inputFile" --L0HTOS --L0GTIS --HLT1
		sbatch SlurmJOBsubmit.sh $inputFile $spl $pol
	done
done
