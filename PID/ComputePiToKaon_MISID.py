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
pathToOutput = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/HistoPID/Pi2K/'
polarities = ['MagUp','MagDown']
stripping = 'Turbo16'
particleList = ['Pi']

Kid = 'IsMuon==0 && DLLK>4.0'
cuts = {'Pi':Kid}


#\"[DLLK > 0.0, DLLK > 4.0 && DLLp < 0.0]\" 
for polarity in polarities:
    for particle in particleList:
        command = 'python '+pathToScript+' \"'+ stripping+'\" \"'+ polarity +'\"  \"'+ particle+ '\" \"[ ' + cuts[particle] + ']\" \"P\" \"ETA\" \"nTracks_Brunel\" -o '+pathToOutput
        print command
        subprocess.call(command,shell=True)


#python $PIDPERFSCRIPTSROOT/scripts/python/MultiTrack/MakePerfHistsRunRange.py  '20' "MagUp" "K" "[DLLK > 0.0]" "P" "ETA" "nTracks"
