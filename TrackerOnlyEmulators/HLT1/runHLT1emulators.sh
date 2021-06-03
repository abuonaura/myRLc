#!/bin/bash

samples=("Lb_Lcmunu" "Lb_Lctaunu" "Lb_LcDs" "Lb_Lc2625munu" "Lb_Lc2625taunu" "Lb_Lc2625Ds" "Lb_Lc2593munu" "Lb_Lc2593taunu" "Lb_Lc2593Ds")
polarities=("MagUp" "MagDown")

for sample in "${samples[@]}";
do
	echo "$sample"
	for polarity in "${polarities[@]}";
	do
		echo "$polarity"
		inputFile="/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/"$sample"_"$polarity".root"
		echo "$inputFile"
		files=()
		for i in {0..10};
		do
			echo "$i"
			start0=$(($i*2000000))
			stop0=$((($i+1)*2000000))
			python emulateHLT1_twoTrack_RDataFrame.py $inputFile $start0 $stop0 $i
			outfile="/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/"$sample"_"$polarity"_wHLT1TwoTracksEmulation_"$i".root"
			echo "$outfile"
			if test -f $outfile; then
				files+=("$outfile")
			fi
		done

		echo ${files[*]}
		finalfile="/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/"$sample"_"$polarity"_wHLT1TwoTracksEmulation.root"
		echo "$finalfile"
		hadd -f $finalfile ${files[*]}
	done
done
