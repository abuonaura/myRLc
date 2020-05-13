"""
Author: Annarita Buonaura
Date: March 1st, 2018

Description:
Compute the probability of wrongly identifying a hadron as a muon => P(h|mu)

How to run it:
1. Make sure to run when starting the session:
     lb-run -c best LHCbDirac/prod bash
     lhcb-proxy-init
2. cd UraniaDev_v7r0
3. ./run bash --norc


Look for esclamation marks for things to fix.
#4.03.20 -> Added DLLe<1 to reduce electron misid
"""
import subprocess

pathToScript = '~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py'
pathToOutput = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/HistoPID/MISIDeff/'
polarities = ['MagUp','MagDown']
stripping = 'Turbo16'
#stripping = '20'
particleList = ['K','Pi','P','Mu']
#particleList = ['Pi']
#cuts = '[DLLmu>-200 && (DLLK<2.0 || DLLK >4.0) && DLLp>0.0 && (MC15TuneV1_ProbNNp - MC15TuneV1_ProbNNK)>0. && IsMuon!=1.0]'

muid = 'IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0 && DLLe<1'
#muid = 'IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0'
cuts = {'K': muid,'Pi':muid,'P': muid, 'Mu':muid}


for polarity in polarities:
    for particle in particleList:
        command = 'python '+pathToScript+' \"'+ stripping+'\" \"'+ polarity +'\"  \"'+ particle+ '\" \"[ ' + cuts[particle] + ']\" \"P\" \"ETA\" \"nTracks_Brunel\" -o '+pathToOutput
        print command
        subprocess.call(command,shell=True)


#python $PIDPERFSCRIPTSROOT/scripts/python/MultiTrack/MakePerfHistsRunRange.py  '20' "MagUp" "K" "[DLLK > 0.0]" "P" "ETA" "nTracks"
