'''
author: A. Buonaura
Date: 15th January 2019
Comment: This scripts compare the ISOLATION BDT distributions for Lb_LcDs and Lb_Lctaunu and try to define a cut to separate signal and bkg region.
'''

import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas, TMath, TLegend

from array import array

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MC/'
datafiles = {'Signal':datadir+'Lb_Lctaunu_PID_reduced_preselected.root','DoubleCh':datadir+'Lb_LcDs_PID_reduced_preselected.root'}

h_ISOBDT_sig = TH1F('h_ISOBDT_sig','',50,-1,1)
h_ISOBDT_2ch = TH1F('h_ISOBDT_2ch','',50,-1,1)
h_ISOBDT_2chK = TH1F('h_ISOBDT_2chK','',50,-1,1)

f_sig = TFile(datafiles['Signal'],'READ')
t_sig = f_sig.Get('DecayTree')

f_2ch = TFile(datafiles['DoubleCh'],'READ')
t_2ch = f_2ch.Get('DecayTree')

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

for i in range(t_sig.GetEntries()):
    t_sig.GetEntry(i)
    if t_sig.Lb_BKGCAT<50 and t_sig.Lc_BKGCAT<30:
        h_ISOBDT_sig.Fill(t_sig.Lb_ISOLATION_BDT)

for j in range(t_2ch.GetEntries()):
    t_2ch.GetEntry(j)
    if t_2ch.Lb_BKGCAT<50 and t_2ch.Lc_BKGCAT<30:
        h_ISOBDT_2ch.Fill(t_2ch.Lb_ISOLATION_BDT)
        if r.TMath.Abs(t_2ch.Lb_TrueHadron_D2_ID)==321:
            h_ISOBDT_2chK.Fill(t_2ch.Lb_ISOLATION_BDT)
                
h_ISOBDT_sig = ScaleHisto(h_ISOBDT_sig,1)
h_ISOBDT_2ch = ScaleHisto(h_ISOBDT_2ch,1)
h_ISOBDT_2chK = ScaleHisto(h_ISOBDT_2chK,1)

legend = r.TLegend(0.6,0.75,0.9,0.9)
legend.AddEntry(h_ISOBDT_sig,'Signal','l')
legend.AddEntry(h_ISOBDT_2ch,'Double Charm','l')

r.gStyle.SetOptStat(0)
c = TCanvas('c','',500,500)
h_ISOBDT_sig.Draw('histo')
h_ISOBDT_sig.GetYaxis().SetRangeUser(0,0.12)
h_ISOBDT_sig.GetXaxis().SetTitle("ISOLATION_BDT")
h_ISOBDT_2ch.SetLineColor(r.kRed)
h_ISOBDT_2ch.Draw('histo sames')
legend.Draw()
c.SaveAs('plots/ISOBDTcomparison.png')

legend0 = r.TLegend(0.6,0.75,0.9,0.9)
legend0.AddEntry(h_ISOBDT_sig,'Signal','l')
legend0.AddEntry(h_ISOBDT_2ch,'Double Charm','l')
legend0.AddEntry(h_ISOBDT_2chK,'Double Charm with K id cut','l')

c0 = TCanvas('c0','',500,500)
h_ISOBDT_sig.Draw('histo')
h_ISOBDT_sig.GetYaxis().SetRangeUser(0,0.2)
h_ISOBDT_2ch.SetLineColor(r.kRed)
h_ISOBDT_2ch.Draw('histo sames')
h_ISOBDT_2chK.SetLineColor(r.kRed)
h_ISOBDT_2chK.SetLineStyle(2)
h_ISOBDT_2chK.Draw('histo sames')
legend0.Draw()
c0.SaveAs('plots/ISOBDTcomparisonWithK.png')


ISOBDT = r.RooRealVar("ISOBDT",'ISOBDT',-1.,1)
sig = r.RooDataHist('sig','ISOBDT Lb->LcTauNu',r.RooArgList(ISOBDT),h_ISOBDT_sig)
charm = r.RooDataHist('charm','ISOBDT Lb->LcXc',r.RooArgList(ISOBDT),h_ISOBDT_2ch)
n_sig = sig.sumEntries()
n_charm = charm.sumEntries()

isobdtcut = []
eff_sig = []
eff_charm = []
significance = []

ndiv = 20

for i in range(ndiv):
    isobdtcut.append(-0.20+i*0.04)
    cut= "ISOBDT<"+str(isobdtcut[i])
    nsel_sig = sig.sumEntries(cut)
    nsel_charm = charm.sumEntries(cut)
    eff_sig.append(nsel_sig*1./n_sig)
    eff_charm.append(nsel_charm*1./n_charm)
    significance.append(nsel_sig/TMath.Sqrt(nsel_sig+nsel_charm))


xgraph = array("d", isobdtcut)
ygraph = array("d",eff_sig)
ygraph1 = array("d",eff_charm)


g_sig = r.TGraph(ndiv,xgraph,ygraph)
g_charm = r.TGraph(ndiv,xgraph,ygraph1)
c1 = TCanvas('c1','c1',500,500)
g_sig.SetTitle('Events with ISOBDT<cut')
g_sig.GetXaxis().SetTitle('efficiency')
g_sig.GetYaxis().SetRangeUser(0,1)
g_sig.SetMarkerColor(r.kRed)
g_sig.Draw('A*')
g_charm.Draw('*')
legend1=r.TLegend(0.6,0.75,0.9,0.9)
legend1.AddEntry(g_sig,'Selected evts signal', "p")
legend1.AddEntry(g_charm,'Selected evts charm', "p")
legend1.Draw()
c1.SaveAs('plots/SelectionEfficiency.png')

ygraph2 = array("d",significance)
g_significance = r.TGraph(ndiv,xgraph,ygraph2)
g_significance.SetTitle('S/#Sqrt(S+B)')
c2 = TCanvas('c2','c2',500,500)
g_significance.Draw('A*')
c2.SaveAs('plots/Significance.png')

