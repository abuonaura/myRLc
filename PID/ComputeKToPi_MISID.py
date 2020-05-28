"""
Author: Annarita Buonaura
Date: March 1st, 2018

Description:
Compute PID efficiencies

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
pathToOutput = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/HistoPID/K2Pi/'
polarities = ['MagUp','MagDown']
stripping = 'Turbo16'
#stripping = '20'
#particleList = ['K','Pi','P','Mu']
particleList = ['K']
#cuts = '[DLLmu>-200 && (DLLK<2.0 || DLLK >4.0) && DLLp>0.0 && (MC15TuneV1_ProbNNp - MC15TuneV1_ProbNNK)>0. && IsMuon!=1.0]'


cuts = {'K':'IsMuon==0 && DLLK<2.0'}


#\"[DLLK > 0.0, DLLK > 4.0 && DLLp < 0.0]\" 
for polarity in polarities:
    for particle in particleList:
        filename = 'Lb_Data_FakeMu_'+polarity+'_first5files_eff_'+particle+'.root'
        command = 'python '+pathToScript+' \"'+ stripping+'\" \"'+ polarity +'\"  \"'+ particle+ '\" \"[ ' + cuts[particle] + ']\" \"P\" \"ETA\" \"nTracks_Brunel\" -o '+pathToOutput
        print command
        subprocess.call(command,shell=True)


#python $PIDPERFSCRIPTSROOT/scripts/python/MultiTrack/MakePerfHistsRunRange.py  '20' "MagUp" "K" "[DLLK > 0.0]" "P" "ETA" "nTracks"