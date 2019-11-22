import ROOT as r
import sys, os
import argparse
from ReduceTree import *
from TrainBDT import *
from CreateFinalWorkingSamples import *
from AddSweights import *
#from IsolationSplit import *

parser = argparse.ArgumentParser(description='Run Preselection on data and MC samples')
parser.add_argument('datatype',choices=['Lctaunu','Lcmunu','LcDs','Lc2593munu','Lc2593taunu','Lc2625munu','Lc2625taunu','all'], help = 'which mc sample we want to run on')
parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on')
parser.add_argument('step',choices=['Trigger&Mcut','addBDTinfo','preselection','all'], help='which step of the preselection chain we want to run')
args = parser.parse_args()

datatype=args.datatype
polarity=args.polarity
step=args.step

datadir = '$FILEDIR/'

datatypes=['Lctaunu','Lcmunu','LcDs','Lc2593munu','Lc2593taunu','Lc2625munu','Lc2625taunu']
#datatypes=['Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2625munu','Lc2625taunu']
polarities=['MagUp','MagDown']
steps=['Trigger&Mcut','addBDTinfo','preselection']

if datatype!='all':
    datatypes=[datatype]
if polarity!='all':
    polarities=[polarity]
if step!='all':
    steps=[step]

print(datatypes)
print(polarities)
print(steps)

BDTcut=0.65

for dtype in datatypes:
    print(dtype)
    if polarity=='all':
        startfname = datadir + 'MC/Lb_'+dtype+'_PID.root'
    else:
        startfname = datadir + 'MC/Lb_'+dtype+'_'+polarity+'_PID.root'
        for stp in steps:
            print(stp)
            if stp=='Trigger&Mcut':
                ofname = startfname[0:-5]+'_reduced.root'
                print(ofname)
                ReduceTree(startfname, 'DecayTree', ofname, 'DecayTree', 'MC')
            if stp=='addBDTinfo':
                model=RunBDT(BDTcut)
                #fname = datadir + 'MC/Lb_'+dtype+'_PID_reduced.root'
                fname = startfname[0:-5]+'_reduced.root'
                print(fname)
                AddBDTinfo(fname,'MC',model)
            if stp=='preselection':
                #fname = datadir + 'MC/Lb_'+dtype+'_PID_reduced.root'
                fname = startfname[0:-5]+'_reduced.root'
                print(fname)
                ApplyFinalSelections(fname,'MC',BDTcut, dtype)
