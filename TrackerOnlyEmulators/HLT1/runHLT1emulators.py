import sys, os, time, subprocess
import emulateHLT1_twoTrack_RDataFrame

inputFile = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/Lb_Lcmunu_MagUp.root'
for i in range(0,20):
    print(i)
    emulateHLT1_twoTrack_RDataFrame.main(inputFile,i*1000000,(i+1)*1000000,i)
    process = subprocess.Popen(emulateHLT1_twoTrack_RDataFrame)
    process.terminate()
    time.sleep(10)
