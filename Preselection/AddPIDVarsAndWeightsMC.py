'''
NOTE: This code must be run in an Urania environment!!!
'''


import sys, os
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_PIDGenWghts_PIDCalibWghts import AddPIDGenWeights
from Add_PIDGenWghts_PIDCalibWghts import AddPIDCalibWeights

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'
UraniaDir = '/home/hep/buonaura/UraniaDev_v7r0/'

if __name__ == "__main__":
    tname = 'tupleout/DecayTree'
    #tname = 'DecayTree'
    #dtype = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
    #dtype = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625Ds']
    dtype = ['Lc2625taunu']
    #polarities=['MagUp','MagDown']
    polarities=['MagUp']
    for polarity in polarities:
        print(polarity)
        for dt in dtype:
            print('   ', dt)
            inputFile = filedir+'Lb_'+dt+'_'+polarity+'_full.root'
            outFilePIDGen = inputFile[0:-5]+'_PIDGen.root'
            outFilePIDCalib = inputFile[0:-5]+'_PIDCalib.root'
            AddPIDGenWeights(inputFile, tname, outFilePIDGen, polarity,UraniaDir)
            AddPIDCalibWeights(inputFile, tname, outFilePIDCalib, polarity)

