#!/bin/bash

path="/home/hep/amathad/Packages/UraniaDev_v8r0"
pidpath=$path"/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack"
perfhistpath="/disk/lhcb_data/amathad/Lb2Lclnu_analysis/perfhist/RLc" #directory to store performance histogram

mkdir -p log

for stripp in "Turbo16"
do 
        for magtype in "MagUp" "MagDown"
        do  
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Mu \[DLLmu\>2\&\&DLLmu-DLLK\>2\&\&DLLmu-DLLp\>2\&\&IsMuon\=\=1\&\&MC15TuneV1_ProbNNghost\<0.2\] Brunel_P Brunel_PT nTracks_Brunel -c Brunel_P\>3000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>16\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b $perfhistpath/binning-"$stripp".py -s binning-Mu-"$stripp"-"$magtype" -q -T > log/Mu-"$stripp"-"$magtype"-Ghost.log
		
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Pi \[DLLK\<2\] Brunel_P Brunel_PT nTracks_Brunel -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b $perfhistpath/binning-"$stripp".py -s binning-Pi-"$stripp"-"$magtype" -q -T > log/Pi-"$stripp"-"$magtype".log

                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype K \[DLLK\>4\] Brunel_P Brunel_PT nTracks_Brunel -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b $perfhistpath/binning-"$stripp".py -s binning-K-"$stripp"-"$magtype" -q -T > log/K-"$stripp"-"$magtype".log

                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype P \[DLLp\>0\] Brunel_P Brunel_PT nTracks_Brunel -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b $perfhistpath/binning-"$stripp".py -s binning-P-"$stripp"-"$magtype" -q -T > log/P-"$stripp"-"$magtype".log
        done
done

unset path
unset pidpath
#screen -dm bash -c "CommandToRun"
