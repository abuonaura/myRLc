import ROOT as r
import sys, os
import argparse
from ReduceTree import *
from TrainBDT import *
from CreateFinalWorkingSamples import *
from AddSweights import *
#from IsolationSplit import *

parser = argparse.ArgumentParser(description='Run Preselection on data and MC samples')
parser.add_argument('datatype',choices=['Data','FakeMu','FakeMuSS','DataSS','all'], help = 'which data sample we want to run on')
parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on')
parser.add_argument('step',choices=['Trigger&Mcut','addBDTinfo','preselection','ISOsplit','all'], help='which step of the preselection chain we want to run')
parser.add_argument('category', choices=['full','iso','all'], help ='which isolation category we want to process')
args = parser.parse_args()

datatype=args.datatype
polarity=args.polarity
step=args.step
category=args.category

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'

datatypes=['Data','FakeMu','FakeMuSS','DataSS']
polarities=['MagUp','MagDown']
steps=['Trigger&Mcut','addBDTinfo','preselection','ISOsplit']


if datatype!='all':
    datatypes=[datatype]
if polarity!='all':
    polarities=[polarity]
if step!='all':
    steps=[step]
if category=='all':
    categories = ['full','iso']
else:
    categories = [category]
                
BDTcut=0.65

for dtype in datatypes:
    print(dtype)
    for pol in polarities:
        startfname = datadir + 'Data/Lb_'+dtype+'_'+pol+'.root'
        for stp in steps:
            if stp=='Trigger&Mcut':
                ofname = datadir + 'Data/Lb_'+dtype+'_'+pol+'_reduced.root'
                ReduceTree(startfname, 'tupleout/DecayTree', ofname, 'DecayTree', dtype)
            if stp=='addBDTinfo':
                model=RunBDT(BDTcut)
                fname = datadir + 'Data/Lb_'+dtype+'_'+pol+'_reduced.root'
                AddBDTinfo(fname,dtype,model)
            if stp=='preselection':
                fname = datadir + 'Data/Lb_'+dtype+'_'+pol+'_reduced.root'
                ApplyFinalSelections(fname,dtype,BDTcut,'')
            if stp=='ISOsplit':
                for cat in categories:
                    os.system('python3 IsolationSplit.py %s %s %s'%(dtype, pol, category))



'''
polarities = ['MagUp']
#mcsamples = ['LcTau','LcMu','LcDs','Lc2593Mu','Lc2593Tau','Lc2625Mu','Lc2625Tau']
mcsamples = ['Lc2625Tau']



startfile = {'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp.root','MagDown':datadir+'Data/Lb_Data_MagDown.root'},
        'DataSS':{'MagUp':datadir+'Data/Lb_DataSS_MagUp.root','MagDown':datadir+'Data/Lb_DataSS_MagDown.root'},
        'FakeMu':{'MagUp':datadir+'Data/Lb_FakeMu_MagUp.root',
                  'MagDown':datadir+'Data/Lb_FakeMu_MagDown.root'},
        'FakeMuSS':{'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp.root',
                  'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown.root'},
        'MC':{'LcTau':datadir+'MC/Lb_Lctaunu_PID.root',
              'LcMu':datadir+'MC/Lb_Lcmunu_PID.root',
              'LcDs':datadir+'MC/Lb_LcDs_PID.root',
              'Lc2593Mu':datadir+'MC/Lb_Lc2593munu_PID.root',
              'Lc2593Tau':datadir+'MC/Lb_Lc2593taunu_PID.root',
              'Lc2625Mu':datadir+'MC/Lb_Lc2625munu_PID.root',
              'Lc2625Tau':datadir+'MC/Lb_Lc2625taunu_PID.root'}}


reducedfile = {'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced.root','MagDown':datadir+'Data/Lb_Data_MagDown_reduced.root'},
        'DataSS':{'MagUp':datadir+'Data/Lb_DataSS_MagUp_reduced.root','MagDown':datadir+'Data/Lb_DataSS_MagDown_reduced.root'},
        'FakeMu':{'MagUp':datadir+'Data/Lb_FakeMu_MagUp_reduced.root',
                  'MagDown':datadir+'Data/Lb_FakeMu_MagDown_reduced.root'},
        'FakeMuSS':{'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp_reduced.root',
                  'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown_reduced.root'},
        'MC':{'LcTau':datadir+'MC/Lb_Lctaunu_PID_reduced.root',
              'LcMu':datadir+'MC/Lb_Lcmunu_PID_reduced.root',
              'LcDs':datadir+'MC/Lb_LcDs_PID_reduced.root',
              'Lc2593Mu':datadir+'MC/Lb_Lc2593munu_PID_reduced.root',
              'Lc2593Tau':datadir+'MC/Lb_Lc2593taunu_PID_reduced.root',
              'Lc2625Mu':datadir+'MC/Lb_Lc2625munu_PID_reduced.root',
              'Lc2625Tau':datadir+'MC/Lb_Lc2625taunu_PID_reduced.root',
              }}


for dtype in datatype:
    if dtype!='MC':
        for polarity in polarities:
            print()
            print('>>>> Processing ', startfile[dtype][polarity])
            print()
            ReduceTree(startfile[dtype][polarity],'tupleout/DecayTree',reducedfile[dtype][polarity],'DecayTree',dtype)
    else:
        for sample in mcsamples:
            print()
            print('>>>> Processing ', startfile[dtype][sample])
            print()
            ReduceTree(startfile[dtype][sample],'DecayTree',reducedfile[dtype][sample],'DecayTree',dtype)

BDTcut = 0.65
#model = RunBDT(BDTcut)
preselfname = {'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced_preselected.root',
                  'MagDown':datadir+'Data/Lb_Data_MagDown_reduced_preselected.root'},
        'DataSS':{'MagUp':datadir+'Data/Lb_DataSS_MagUp_reduced_preselected.root',
                  'MagDown':datadir+'Data/Lb_DataSS_MagDown_reduced_preselected.root'},
        'FakeMu':{'MagUp':datadir+'Data/Lb_FakeMu_MagUp_reduced_preselected.root',
                  'MagDown':datadir+'Data/Lb_FakeMu_MagDown_reduced_preselected.root'},
        'FakeMuSS':{'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp_reduced_preselected.root',
                  'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown_reduced_preselected.root'}}

for dtype in datatype:
    if dtype!='MC':
        for polarity in polarities:
            print()
            print('>>>> Processing ', reducedfile[dtype][polarity])
            print()
            #AddBDTinfo(reducedfile[dtype][polarity],dtype,model)
            #if dtype=='Data' and polarity == 'MagUp':
                #bdtfname = ApplyBDTcut(reducedfile[dtype][polarity], dtype,BDTcut)
                #swfname = ComputeSweights(bdtfname,dtype,polarity)
                #os.system('python3 ReconstructDpDs.py')
            preselfname = ApplyFinalSelections(reducedfile[dtype][polarity],dtype,BDTcut)
            #CreateSplitSamples(preselfname[dtype][polarity],dtype,polarity)
                
    else:
        for sample in mcsamples:
            print()
            print('>>>> Processing ', reducedfile[dtype][sample])
            print()
            AddBDTinfo(reducedfile[dtype][sample],dtype,model)
            preselfname = ApplyFinalSelections(reducedfile[dtype][sample],dtype,BDTcut)
'''

