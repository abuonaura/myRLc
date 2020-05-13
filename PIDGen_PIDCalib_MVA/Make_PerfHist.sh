#!/bin/bash

path=$1
#path="/home/hep/amathad/Packages/UraniaDev_v8r0"
pidpath=path"PIDCalib/PIDPerfScripts/scripts/python/MultiTrack"

perfhistpath=$2
#perfhistpath="/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist"

mkdir -p log

#run2
for stripp in "Turbo16"
do 
        for magtype in "MagUp" "MagDown"
        do  
                echo "bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Mu \[DLLmu\>2\&\&DLLmu-DLLK\>2\&\&DLLmu-DLLp\>2\&\&IsMuon\=\=1\] Brunel_P Brunel_PT nTracks -c Brunel_P\>3000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>16\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-Mu-"$stripp"-"$magtype" -q -T > log/Mu-"$stripp"-"$magtype".log 2>&1 &"
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Mu \[DLLmu\>2\&\&DLLmu-DLLK\>2\&\&DLLmu-DLLp\>2\&\&IsMuon\=\=1\] Brunel_P Brunel_PT nTracks -c Brunel_P\>3000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>16\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-Mu-"$stripp"-"$magtype" -q -T > log/Mu-"$stripp"-"$magtype".log 2>&1 &
                echo "bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Pi \[DLLK\<2\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-Pi-"$stripp"-"$magtype" -q -T > log/Pi-"$stripp"-"$magtype".log 2>&1 &"
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype Pi \[DLLK\<2\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-Pi-"$stripp"-"$magtype" -q -T > log/Pi-"$stripp"-"$magtype".log 2>&1 &
                echo "bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype K \[DLLK\>4\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-K-"$stripp"-"$magtype" -q -T > log/K-"$stripp"-"$magtype".log 2>&1 &"
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype K \[DLLK\>4\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-K-"$stripp"-"$magtype" -q -T > log/K-"$stripp"-"$magtype".log 2>&1 &
                echo "bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype P \[DLLp\>0\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-P-"$stripp"-"$magtype" -q -T > log/P-"$stripp"-"$magtype".log 2>&1 &"
                bash $path/run python $pidpath/MakePerfHistsRunRange.py $stripp $magtype P \[DLLp\>0\] Brunel_P Brunel_PT nTracks -c Brunel_PT\>300\&\&Brunel_P\>2000\&\&nSPDHits\<600\&\&Brunel_IPCHI2\>9\&\&Brunel_TRACK_GHOSTPROB\<0.5 -o $perfhistpath -b /home/hep/amathad/LbToLclnu_2016/PID_RLc/binoutput/binning-"$stripp".py -s binning-P-"$stripp"-"$magtype" -q -T > log/P-"$stripp"-"$magtype".log 2>&1 &
        done
done

unset path
unset pidpath
#screen -dm bash -c "CommandToRun"
