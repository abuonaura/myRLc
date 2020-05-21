'''
Author: Annarita Buonaura
Date: May 21, 2020

Description: Script to run Abhijit's code to create 2 root files:
    1. the one named with "_PIDGen": p_PROBNNp_corr, K_PROBNNK_corr, pi_PROBNNpi_corr
    2. the one named with "_PIDCalib": PIDCalib weights (name of the variable: Event_PIDCalibEffWeight)


How to run:
    NOTE: This code must be run in an Urania environment!!!

    - On full simulation (all samples, producing from scratch also the trigger files):
        nohup python -i AddPIDVarsAndWeightsMC.py --MCfull all all > AddPID.txt &

    - On TrackerOnly simulation (all samples, producing from scratch also the trigger files):
	 nohup python -i AddPIDVarsAndWeightsMC.py --MCTrackerOnly all all > AddPID.txt &

    - For one sample and one polarity:
        nohup python -i AddPIDVarsAndWeightsMC.py --MCTrackerOnly Lcmunu MagUp > AddPID.txt &

'''


import sys, os
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_PIDGenWghts_PIDCalibWghts import AddPIDGenWeights
from Add_PIDGenWghts_PIDCalibWghts import AddPIDCalibWeights
from argparse import ArgumentParser

UraniaDir = '/home/hep/buonaura/UraniaDev_v7r0/'

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    parser.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    parser.add_argument('datatype',choices=['Lctaunu','Lcmunu','LcDs','Lc2593munu','Lc2593taunu','Lc2625munu','Lc2625taunu','all'], help = 'which mc sample we want to run on', default='all')
    parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on', default = 'all')

    options = parser.parse_args()
    MCfull = options.MCfull
    MCTO = options.MCTO
    datatype=options.datatype
    polarity=options.polarity

    if MCfull==True:
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
    if MCTO==True:
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/'

    tname = 'tupleout/DecayTree'

    datatypes = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
    polarities=['MagUp','MagDown']

    if datatype!='all':
        datatypes=[datatype]
    if polarity!='all':
        polarities=[polarity]

    for polarity in polarities:
        print(polarity)
        for dt in datatypes:
            print('   ', dt)
	    if MCfull==True:
                inputFile = filedir+'Lb_'+dt+'_'+polarity+'_full.root'
            if MCTO==True:
                inputFile = filedir+'Lb_'+dt+'_'+polarity+'.root'
            outFilePIDGen = inputFile[0:-5]+'_PIDGen.root'
            outFilePIDCalib = inputFile[0:-5]+'_PIDCalib.root'
            AddPIDGenWeights(inputFile, tname, outFilePIDGen, polarity,UraniaDir)
            AddPIDCalibWeights(inputFile, tname, outFilePIDCalib, polarity)

