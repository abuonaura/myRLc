'''
author: A.Buonaura
date: 7 January 2018
This scripts plots different Lc_M distributions (standard, with sweights, with muPIDweight, with muPIDweight+sweights).
In the last plot it compares the standard and the muPIDweigthed Lc_M distributions. Being similar, it makes no sense to compute the sweights before or after having applied the muPIDweight.
'''

import ROOT as r
import os.path
import os,sys,getopt,time

datadir = '$FILEDIR/MISID/SameSign/'

suffix = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

if __name__ == '__main__':
    opts, args = getopt.getopt(sys.argv[1:], "",["full","iso","Kenriched"])
    print (opts,args)
    for o, a in opts:
        if o in ("--full",):
            sample = 'full'
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'

    f = r.TFile(datadir+'K_sample_MagDown'+suffix[sample],'READ')
    t = f.Get('DecayTree')

    h=r.TH1F('h','',100,2230,2330)
    h1=r.TH1F('h1','',100,2230,2330)
    h2=r.TH1F('h2','',100,2230,2330)
    h3=r.TH1F('h3','',100,2230,2330)
    t.Draw('Lc_M>>h','','goff')
    t.Draw('Lc_M>>h1','sw_sig*(Lc_M>0)','goff')
    t.Draw('Lc_M>>h2','w_recomu_CF*(Lc_M>0)','goff')
    t.Draw('Lc_M>>h3','sw_sig*w_recomu_CF*(Lc_M>0)','goff')

    c = r.TCanvas('c','c',500,500)
    h1.Draw('histo')
    h1.SetLineColor(r.kRed)
    h.Draw('histo sames')

    c1 = r.TCanvas('c1','c1',500,500)
    h2.SetLineColor(r.kBlack)
    h2.Draw('histo')
    c1a = r.TCanvas('c1a','c1a',500,500)
    h3.SetLineColor(r.kRed)
    h3.Draw('histo')


    '''
    Here you compare the shapes of the Lc_M peak before and after applying the PID weight (=prob of wongly reco the hadron as a muon/prob of correctly identifying the hadron) taking the crossfeed into account. If they are similar, then it means that it is not such a big problem that you computed the sweights at the beginning of the analysis
    '''


    ScaleHisto(h,1)
    ScaleHisto(h2,1)

    c2 = r.TCanvas('c2','c2',500,500)
    h.Draw()
    h2.SetLineColor(r.kRed)
    h2.Draw('sames')
