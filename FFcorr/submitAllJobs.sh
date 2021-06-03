#!/bin/bash

samples=("Lcmunu" "Lctaunu" "Lc2593munu" "Lc2625munu" "Lc2593taunu" "Lc2625taunu")
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
