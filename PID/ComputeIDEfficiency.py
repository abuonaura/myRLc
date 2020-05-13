"""
Author: Annarita Buonaura
Date: March 1st, 2018

Description:
Computes the efficiency of correctly identifying hadrons and muons => P(h|h)

How to run it:
1. Make sure to run when starting the session:
     lb-run -c best LHCbDirac/prod bash
     lhcb-proxy-init
2. cd UraniaDev_v7r0
3. ./run bash --norc


Look for esclamation marks for things to fix.
"""
import subprocess

pathToScript = '~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py'
pathToOutput = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/HistoPID/IDeff/'
polarities = ['MagUp','MagDown']
stripping = 'Turbo16'
#stripping = '20'
particleList = ['K','Pi','P','Mu']
#particleList = ['K']
#cuts = '[DLLmu>-200 && (DLLK<2.0 || DLLK >4.0) && DLLp>0.0 && (MC15TuneV1_ProbNNp - MC15TuneV1_ProbNNK)>0. && IsMuon!=1.0]'
muid = 'IsMuon==0'
cuts = {'K':'DLLK >4.0 &&'+muid,'Pi':'DLLK<2 && '+muid,'P':'DLLp < 0.0 && (MC15TuneV1_ProbNNp - MC15TuneV1_ProbNNK)>0. &&' + muid, 'Mu':'DLLmu>-200 &&'+muid}
#\"[DLLK > 0.0, DLLK > 4.0 && DLLp < 0.0]\" 
for polarity in polarities:
    for particle in particleList:
        filename = 'Lb_Data_FakeMu_'+polarity+'_first5files_eff_'+particle+'.root'
        command = 'python '+pathToScript+' \"'+ stripping+'\" \"'+ polarity +'\"  \"'+ particle+ '\" \"[ ' + cuts[particle] + ']\" \"P\" \"ETA\" \"nTracks_Brunel\"  -o '+pathToOutput
        print command
        subprocess.call(command,shell=True)


#python $PIDPERFSCRIPTSROOT/scripts/python/MultiTrack/MakePerfHistsRunRange.py  '20' "MagUp" "K" "[DLLK > 0.0]" "P" "ETA" "nTracks"
