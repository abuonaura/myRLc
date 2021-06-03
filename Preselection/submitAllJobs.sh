#!/bin/bash

samples=("Lb_Lcmunu" "Lb_Lctaunu" "Lb_LcDs" "Lb_Lc2625munu" "Lb_Lc2625taunu" "Lb_Lc2625Ds" "Lb_Lc2593munu" "Lb_Lc2593taunu" "Lb_Lc2593Ds" "B_Lcpbarmunu" "Lb_Lc2765munu" "Lb_Lc2880munu")
polarities=("MagUp" "MagDown")

for sample in "${samples[@]}";
do
	echo "$sample"
	for polarity in "${polarities[@]}";
	do
		echo "$polarity"
		A=$sample
		B=$polarity
		sbatch --job-name=$A.$B.run --output=logs/$A.$B.out --export=A=$A,B=$B SlurmJOBsubmit.sh $sample $polarity
	done
done
