import ROOT as r
from TruthMatch import *

import argparse

parser = argparse.ArgumentParser(description="Make my DaVinci job.")
parser.add_argument('-s','--sample', choices=['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds','all'],
                                    help='which sample the script is run on', required=False)
parser.add_argument('-p','--polarity', choices=['MagUp', 'MagDown','all'],
                                    help='Polarity of data-taking to run over', required=False)
parser.add_argument('-t','--type', choices=['MC', 'MCtrackeronly'],
                                    help='Polarity of data-taking to run over', required=False)

args = parser.parse_args()

sample = args.sample
polarity = args.polarity
dtype = args.type

if sample!='all':
    samples = [sample]
else:
    samples = ['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds']

if polarity!='all':
    polarities = [polarity]
else:
    polarities = ['MagUp','MagDown']


for s in samples:
    cut = GetTruthMatchCut(s)
    for p in polarities:
        print ('Running on : ', s, p)
        if dtype=='MCtrackeronly':
            fname = '$FILEDIR/MCtrackeronly/'+s+'_'+p+'.root'
            tname = 'tupleout/DecayTree'
        if dtype=='MC':
            fname = '$FILEDIR/MC/'+s+'_'+p+'_PID.root'
            tname = 'DecayTree'

        f = r.TFile(fname,'READ')
        t = f.Get(tname)

        ofname = fname[0:-5]+'_TMatched.root'
        of = r.TFile(ofname,'RECREATE')
        ot = r.TTree('DecayTree','DecayTree')

        ot = t.CopyTree(cut.GetTitle())
        of.Write()
        of.Close()
        f.Close()
