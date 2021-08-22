"""
Author: Annarita Buonaura
Date: August 20, 2021
"""
import subprocess

pathToScript = '~/UraniaDev_v7r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/MakePerfHistsRunRange.py'
pathToOutput = '/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/'
polarities = ['MagUp','MagDown']
stripping = 'Turbo16'
particleList = ['Pi','K']

Kid = 'IsMuon==0 && DLLK>4.0, IsMuon==0 && DLLK>4.0 && DLLp-DLLK<0, IsMuon==0 && DLLK>4.0 && DLLp-DLLK>=0'


for polarity in polarities:
    for particle in particleList:
        command = 'python '+pathToScript+' \"'+ stripping+'\" \"'+ polarity +'\"  \"'+ particle+ '\" \"[ ' + Kid + ']\" \"Brunel_P\" \"Brunel_ETA\" \"nTracks_Brunel\" -o '+pathToOutput
        print(command)
        subprocess.call(command,shell=True)


#python $PIDPERFSCRIPTSROOT/scripts/python/MultiTrack/MakePerfHistsRunRange.py  '20' "MagUp" "K" "[DLLK > 0.0]" "P" "ETA" "nTracks"
