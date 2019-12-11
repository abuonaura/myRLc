import ROOT as r
import argparse
from ReduceTree import *
from CreateFinalWorkingSamples import *
from TruthMatch import *

parser = argparse.ArgumentParser(description='Run Preselection on data and MC samples')
parser.add_argument('sample',choices=['Lctaunu','Lcmunu','LcDs','Lc2593munu','Lc2593taunu','Lc2625munu','Lc2625taunu','Lc2625Ds','Lc2593Ds'], help = 'which mc sample we want to run on')
parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on')
args = parser.parse_args()

sample=args.sample
polarity=args.polarity

if polarity!='all':
    polarities = [polarity]
else:
    polarities = ['MagUp','MagDown']

nstart=0
n_tm=0
n_trigger=0
n_stripPID=0
n_LcM=0
n_BDT=0

for p in polarities:
    f = r.TFile('$FILEDIR/MC/Lb_'+sample+'_'+p+'_PID.root','READ')
    t = f.Get('DecayTree')
    nstart += t.GetEntries()
    cut_tm = GetTruthMatchCut('Lb_'+sample)
    cut_trigger = GetTriggerCut('MC')
    cut_stripPID = GetStrippingPIDcut()
    cut_LcM = r.TCut("Lc_M>2230 && Lc_M<2330")

    n_tm += t.Draw('Lc_M',cut_tm,'goff')
    n_trigger += t.Draw('Lc_M',cut_tm+cut_trigger,'goff')
    n_stripPID += t.Draw('Lc_M',cut_tm+cut_trigger+cut_stripPID,'goff')
    n_LcM += t.Draw('Lc_M',cut_tm+cut_trigger+cut_stripPID+cut_LcM,'goff')

    f1 = r.TFile('$FILEDIR/MC/Lb_'+sample+'_'+p+'_PID_reduced.root','READ')
    t1 = f1.Get('DecayTree')

    cut_BDT = GetBDTcut(0.65)
    n_BDT += t1.Draw('Lc_M',cut_tm+cut_BDT,'goff')

print('Start: ', nstart)
print('After truth matching: {}  eff: {:.3f}'.format(n_tm, n_tm*1./nstart))
print('After trigger: {}  eff: {:.3f}'.format(n_trigger, n_trigger*1./n_tm))
print('After Stripping PID cut: {}  eff: {:.3f}'.format(n_stripPID, n_stripPID*1./n_trigger))
print('After Lc mass cut: {}   eff: {:.3f}'.format(n_LcM, n_LcM*1./n_stripPID))
print('After BDT cut: {}   eff: {:.3f}'.format(n_BDT, n_BDT*1./n_LcM))
