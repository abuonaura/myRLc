#!/bin/python

import os
#import glob
#files=glob.glob(Dir+"*_TMatched.root")

#Evaluate MC
Dir = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/RLc_MC/'
files = [Dir+'Lb_Lc2593taunu_MagUp_TMatched_PID.root']

from Add_PIDGenWghts_PIDCalibWghts import AddPIDGenWeights, AddPIDCalibWeights
UraniaDir = '/home/hep/amathad/Packages/UraniaDev_v8r0/'
for fn in files:
    fn_pidgen   = fn.replace('.root','_PIDGen.root')
    fn_pidcalib = fn.replace('.root','_PIDCalib.root')
    fn_mva      = fn.replace('.root','_MVA.root')
    if 'MagUp' in fn: 
        magtype = 'MagUp'
    else:
        magtype = 'MagDown'

    ##Read ONLY 2000 candidates (for testing purpose only)
    #AddPIDGenWeights(  fn, 'DecayTree', fn_pidgen,   magtype, UraniaDir, nentries_to_read=2000)
    #AddPIDCalibWeights(fn, 'DecayTree', fn_pidcalib, magtype, UraniaDir, nentries_to_read=2000)
    #AddPIDCalibWeights(fn, 'DecayTree', fn_pidcalib, magtype, nentries_to_read=2000)
    AddPIDCalibWeights(fn, 'DecayTree', 'Test.root', magtype, nentries_to_read=2000)

##Make sure all the python packages are present (see README.md)
#rom Add_MVA import AddBDTinfo
#for fn in files:
#    fn_mva      = fn.replace('.root','_MVA.root')
#    AddBDTinfo(fn, 'DecayTree', fn_mva, 'MC', nentries_to_read=2000)
#
##Evaluate Data
#DirData     = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/Data/'
#infnamedata = DirData+'Data_MagDown_2016.root'
#oufnamedata = infnamedata.replace('.root','_MVA.root')
#AddBDTinfo(infnamedata,'tupleout/DecayTree',oufnamedata,'Data', nentries_to_read=2000)
