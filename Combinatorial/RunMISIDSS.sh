#!/bin/sh

output1="subsamples_$1.txt"
output2="crossfeed_$1.txt"
output3="MISIDfractions_$1.txt"
output4="weightMISID_$1.txt"
python3 CreateMISIDSubSamples.py --$1 >> $output1 && echo "Created sub samples" .
if [ -z "$2" ]
then
	python3 EvaluateCrossFeed.py --$1 >> $output2 && echo "Evaluated Cross-feeed" .
else
	python3 EvaluateCrossFeed.py --$1 --$2 >> $output2 && echo "Evaluated Cross-feeed data" .
fi
python3 ComputeMISIDfractions.py --$1 >> $output3 && echo "Computed MISID fractions" .
if [ -z "$2" ]
then
	python3 AddWeightMISIDsubtraction.py --$1 >> $output4 && echo "Added MISID weights" .
else
	python3 AddWeightMISIDsubtraction.py --$1 --$2 >> $output4 && echo "Added MISID weights" .
fi
