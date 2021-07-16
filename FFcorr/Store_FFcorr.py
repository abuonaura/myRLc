#!/bin/python

##Please install following packages before running this
#pip install git+https://github.com/scikit-hep/uproot4.git@master #(uproot4)
#pip install uproot #(uproot3)

### To run:
### - To run on all samples and all polarities, just choose which type of simulation: MCfull/MCTrackerOnly:
###    python3 -i Store_FFcorr.py --MCfull (or -- MCTrackerOnly)
### - To run on one sample and both polarities:
###    python3 -i Store_FFcorr.py --MCfull (or -- MCTrackerOnly) samplename all
### - To run on one sample and one polarity (e.g. MagUp):
###    python3 -i Store_FFcorr.py --MCfull (or -- MCTrackerOnly) samplename MagUp


from Add_FFcorr import AddFFcorr 
import os, sys
import argparse

def init():
    ap = argparse.ArgumentParser(description='Create file with FF correction weights')
    ap.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    ap.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    ap.add_argument('sample', help = 'which mc sample we want to run on', default = 'all')
    ap.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on', default = 'all')
    args = ap.parse_args()
    return args


if __name__== "__main__":
    args = init()
    if args.MCfull==True:
        print('Running over full MC samples')
        MCtype    = 'MCfull'
        directory = '/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/' 
    if args.MCTO==True:
        MCtype = 'MCTrackerOnly'
        directory = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/' 
        print('Running over TrackerOnly MC samples')

    samples    = ['Lcmunu','Lctaunu','Lc2593munu','Lc2625munu','Lc2593taunu','Lc2625taunu']
    polarities = ['MagUp','MagDown']
    Lcstate    = {'Lcmunu':'Lc','Lctaunu':'Lc','Lc2593munu':'Lc2595','Lc2593taunu':'Lc2595','Lc2625munu':'Lc2625','Lc2625taunu':'Lc2625'}
    leptname   = {'Lcmunu':'mu','Lctaunu':'tau','Lc2593munu':'mu','Lc2593taunu':'tau','Lc2625munu':'mu','Lc2625taunu':'tau'}

    if args.sample!='all':
        samples=[args.sample]
    if args.polarity!='all':
        polarities=[args.polarity]

    inputtreename         = 'tupleout/DecayTree'
    outputtreename        = 'DecayTree'
    q2True_branchname     = 'Lb_True_Q2_mu'     #This has to be in MeV^2
    costhlTrue_branchname = 'Lb_True_Costhetal_mu'

    for sample in samples:
        print('Analysing sample: ',sample)
        for polarity in polarities:
            print('Running over polarity: ',polarity)
            if MCtype == 'MCfull':
                fn = directory+'Lb_'+sample+'_'+polarity+'_full.root'
            if MCtype == 'MCTrackerOnly':
                fn = directory+'Lb_'+sample+'_'+polarity+'.root'
            fn_new  = fn[0:-5]+'_FFcorrections.root'
            print(fn, inputtreename, fn_new, outputtreename, leptname, q2True_branchname, costhlTrue_branchname)
            AddFFcorr(fn, inputtreename, fn_new, outputtreename, Lcstate[sample], leptname[sample], q2True_branchname, costhlTrue_branchname, nentries_to_read = 50000000) #only read 1000 cands
    


