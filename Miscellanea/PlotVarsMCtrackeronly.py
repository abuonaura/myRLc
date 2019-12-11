import ROOT as r
import sys
sys.path.append('../Preselection/')
from TruthMatch import *

import argparse

parser = argparse.ArgumentParser(description="Make my DaVinci job.")
parser.add_argument('-s','--sample', choices=['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds','all'],
                                    help='which sample the script is run on', required=False)
parser.add_argument('-p','--polarity', choices=['MagUp', 'MagDown','all'],
                                    help='Polarity of data-taking to run over', required=False)

args = parser.parse_args()

sample = args.sample
polarity = args.polarity

if sample!='all':
    samples = [sample]
else:
    samples = ['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds']

if polarity!='all':
    polarities = [polarity]
else:
    polarities = ['MagUp','MagDown']


colors = {'Lb_Lcmunu':r.kBlue,
          'Lb_Lctaunu': r.kRed,
          'Lb_LcDs': r.kGreen+1,
          'Lb_Lc2625taunu':r.kRed+2,
          'Lb_Lc2625munu':r.kBlue+3,
          'Lb_Lc2593taunu':r.kRed-5,
          'Lb_Lc2593munu':r.kBlue-10,
          'Lb_Lc2625Ds':r.kGreen+3,
          'Lb_Lc2593Ds':r.kGreen-8
        }

nentries={}
tm_entries = {}
tm_evts = {}
eff_tm = {}

def make3Dhisto(sample, polarity,tree,cut):
    hKi3 = r.TH3D("h_"+sample+"_"+polarity,"h_"+sample+"_"+polarity,4,-2,14,10,0,2600,10,-2,14)
    a = tree.Draw("FitVar_Mmiss2_mLc/1.E6:FitVar_El_mLc:FitVar_q2_mLc/1.E6>>h_"+sample+"_"+polarity,cut)
    hKi3 = r.gPad.GetPrimitive("h_"+sample+"_"+polarity)
    return hKi3


histos = {}
for s in samples:
    nentries[s] = {}
    tm_entries[s] = {}
    tm_evts[s] = {}
    eff_tm[s] = {}
    histos[s] = {}
    cut = GetTruthMatchCut(s)
    for p in polarities:
        print ('Running on : ', s, p)
        f = r.TFile('$FILEDIR/MCtrackeronly/'+s+'_'+p+'.root','READ')
        t = f.Get('tupleout/DecayTree')
        nentries[s][p] = t.GetEntries()
        tm_entries[s][p] = t.Draw('Lc_M',cut)
        eff_tm[s][p] = tm_entries[s][p]*1./nentries[s][p]
        histos[s][p]=make3Dhisto(s,p,t,cut)
        histos[s][p].SetDirectory(0)


c = r.TCanvas('c','c',1500,500)
c.Divide(3,1)
#print(len(samples),len(polarities))
#c.Divide(len(samples),len(polarities))
#n=0
hsx = r.THStack("hsx","Stacked q2");
hsy = r.THStack("hsy","Stacked El");
hsz = r.THStack("hsz","Stacked Mmiss2");

l = r.TLegend(0.6,0.5,0.9,0.9)

for s in samples:
    c.cd(1)
    hx = histos[s]['MagUp'].ProjectionX()
    hx.SetLineColor(colors[s])
    hx.SetFillColor(colors[s])
    hx.SetFillStyle(3004)
    hsx.Add(hx)
    hx.Draw('same')
    c.cd(2)
    hy = histos[s]['MagUp'].ProjectionY()
    hy.SetLineColor(colors[s])
    hy.SetFillColor(colors[s])
    hy.SetFillStyle(3004)
    hy.Draw('same')
    hsy.Add(hy)
    c.cd(3)
    hz = histos[s]['MagUp'].ProjectionZ()
    hz.SetLineColor(colors[s])
    hz.SetFillColor(colors[s])
    hz.SetFillStyle(3004)
    hz.Draw('same')
    hsz.Add(hz)
    l.AddEntry(hz,s,'f')
    l.Draw()

r.gStyle.SetOptStat(0)
c1 = r.TCanvas('c1','c1',1500,500)
c1.Divide(3,1)
c1.cd(1)
hsx.Draw()
c1.cd(2)
hsy.Draw()
c1.cd(3)
hsz.Draw()
l.Draw()

#---- print statistics ----

entries = {}
entries_tm = {}
eff = {}
for s in samples:
    entries[s]=0
    entries_tm[s]=0
    for p in polarities:
        entries[s]+=nentries[s][p]
        entries_tm[s]+=tm_entries[s][p]
    eff[s] = entries_tm[s]*1./entries[s]
    print('sample: {} | nentries: {} | after TM: {} | sel_eff: {:.3f}'.format(s,entries[s],entries_tm[s],eff[s]))
