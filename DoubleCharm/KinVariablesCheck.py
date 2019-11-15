'''
author: A. Buonaura
date: 14 January 2019
comment: This routine is meant to familiarise with the LcDs MC file and plot the kinematical distributions Lb_LcDs (2body) and Lb_LcD (3 body)

'''

import ROOT as r
from ROOT import TFile, TTree, TCanvas, TH1F, TLegend
import numpy as np

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MC/'

f = TFile(datadir+'Lb_LcDs_PID_reduced_preselected.root','READ')
t = f.Get('DecayTree')

charm=['D0','D1','D2']

twobody=0
threebody=0

h_E_2b = TH1F('h_E_2b',';E^{*}_{#mu} (GeV);',20,0,2600)
h_E_3b = TH1F('h_E_3b',';E^{*}_{#mu} (GeV);',20,0,2600)
h_q2_2b = TH1F('h_q2_2b',';q^{2} (GeV^{2}/c^{2});;',20,-2E6,18E6)
h_q2_3b = TH1F('h_q2_3b',';q^{2} (GeV^{2}/c^{2});',20,-2E6,18E6)
h_Mmiss2_2b = TH1F('h_Mmiss2_2b',';m^{2}_{miss} (GeV^{2}/c^{4});;',20,-2E6,18E6)
h_Mmiss2_3b = TH1F('h_Mmiss2_3b',';m^{2}_{miss} (GeV^{2}/c^{4});',20,-2E6,18E6)

h_E_0 = TH1F('h_E_0',';E^{*}_{#mu} (GeV);',20,0,2600)
h_q2_0 = TH1F('h_q2_0',';q^{2} (GeV^{2}/c^{2});;',20,-2E6,18E6)
h_Mmiss2_0 = TH1F('h_Mmiss2_0',';q^{2} (GeV^{2}/c^{2});;',20,-2E6,18E6)

h_D0_ID = TH1F('h_D0_ID','PID #Lambda_{b} D0',10000,-5000,5000)
h_D1_ID_1 = TH1F('h_D1_ID_1','PID #Lambda_{b} D1',18,1,19)
h_D1_ID_2 = TH1F('h_D1_ID_2','PID #Lambda_{b} D1',18,1,19)
D1particles = ['D^{+}','D^{0}','D^{*}_{0}(2400)^{+}','D^{*}_{0}(2400)^{0}','D^{*}(2010)^{+}','D^{*}(2007)^{0}','D_{1}(2420)^{+}','D_{1}(2420)^{0}','D_{1}(H)^{+}','D_{1}(2430)^{0}','D^{*}_{2}(2460)^{+}','D^{*}_{2}(2460)^{0}','D^{+}_{s}','D^{*}_{s0}(2317)^{+}','D^{*+}_{s}','D_{s1}(2356)^{+}','D_{s1}(2460)^{+}','D^{*}_{s}(2573)^{+}']
for n, d1p in enumerate(D1particles):
    print(n, d1p)
    h_D1_ID_1.GetXaxis().SetBinLabel(n+1,d1p)
    h_D1_ID_2.GetXaxis().SetBinLabel(n+1,d1p)
h_D1_ID_1.GetXaxis().LabelsOption('v')
h_D1_ID_1.GetXaxis().SetLabelSize(0.05)
h_D1_ID_2.GetXaxis().LabelsOption('v')
h_D1_ID_2.GetXaxis().SetLabelSize(0.05)

h_D2_ID = TH1F('h_D2_ID','PID #Lambda_{b} D2',10000,-5000,5000)

def FillHistoD1(h,pdg):
    if r.TMath.Abs(pdg)==411:
        h.Fill(1)
    if r.TMath.Abs(pdg)==421:
        h.Fill(2)
    if r.TMath.Abs(pdg)==10411:
        h.Fill(3)
    if r.TMath.Abs(pdg)==10421:
        h.Fill(4)
    if r.TMath.Abs(pdg)==413:
        h.Fill(5)
    if r.TMath.Abs(pdg)==423:
        h.Fill(6)
    if r.TMath.Abs(pdg)==10413:
        h.Fill(7)
    if r.TMath.Abs(pdg)==10423:
        h.Fill(8)
    if r.TMath.Abs(pdg)==20413:
        h.Fill(9)
    if r.TMath.Abs(pdg)==20423:
        h.Fill(10)
    if r.TMath.Abs(pdg)==415:
        h.Fill(11)
    if r.TMath.Abs(pdg)==425:
        h.Fill(12)
    if r.TMath.Abs(pdg)==431:
        h.Fill(13)
    if r.TMath.Abs(pdg)==10431:
        h.Fill(14)
    if r.TMath.Abs(pdg)==433:
        h.Fill(15)
    if r.TMath.Abs(pdg)==10433:
        h.Fill(16)
    if r.TMath.Abs(pdg)==20433:
        h.Fill(17)
    if r.TMath.Abs(pdg)==435:
        h.Fill(18)
    h.SetDirectory(0)
    return h

for i in range(t.GetEntries()):
#for i in range(100):
    t.GetEntry(i)
    '''
    print('Lb_TrueHadron_D0_ID = ',t.Lb_TrueHadron_D0_ID)
    print('Lb_TrueHadron_D1_ID = ',t.Lb_TrueHadron_D1_ID)
    print('Lb_TrueHadron_D2_ID = ',t.Lb_TrueHadron_D2_ID)
    print('Lb_TrueHadron_D0_GD0_ID = ',t.Lb_TrueHadron_D0_GD0_ID)
    print('Lb_TrueHadron_D0_GD1_ID = ',t.Lb_TrueHadron_D0_GD1_ID)
    print('Lb_TrueHadron_D0_GD2_ID = ',t.Lb_TrueHadron_D0_GD2_ID)
    print('Lb_TrueHadron_D1_GD0_ID = ',t.Lb_TrueHadron_D1_GD0_ID)
    print('Lb_TrueHadron_D1_GD1_ID = ',t.Lb_TrueHadron_D1_GD1_ID)
    print('Lb_TrueHadron_D1_GD2_ID = ',t.Lb_TrueHadron_D1_GD2_ID)
    print('Lb_TrueHadron_D2_GD0_ID = ',t.Lb_TrueHadron_D2_GD0_ID)
    print('Lb_TrueHadron_D2_GD1_ID = ',t.Lb_TrueHadron_D2_GD1_ID)
    print('Lb_TrueHadron_D2_GD2_ID = ',t.Lb_TrueHadron_D2_GD2_ID)
    print('')
    '''
    if t.Lb_TrueHadron_D0_ID!=0 and t.Lb_TrueHadron_D1_ID!=0 and t.Lb_TrueHadron_D2_ID==0:
        twobody+=1
        h_E_2b.Fill(t.FitVar_El_mLc)
        h_q2_2b.Fill(t.FitVar_q2_mLc)
        h_Mmiss2_2b.Fill(t.FitVar_Mmiss2_mLc)
        h_D0_ID.Fill(t.Lb_TrueHadron_D0_ID)
        h_D1_1=FillHistoD1(h_D1_ID_1,t.Lb_TrueHadron_D1_ID)
    if t.Lb_TrueHadron_D0_ID!=0 and t.Lb_TrueHadron_D1_ID!=0 and t.Lb_TrueHadron_D2_ID!=0:
        threebody+=1
        h_E_3b.Fill(t.FitVar_El_mLc)
        h_q2_3b.Fill(t.FitVar_q2_mLc)
        #h_Mmiss2_2b.Fill(t.FitVar_Mmiss2_mLc)
        h_Mmiss2_3b.Fill(t.FitVar_Mmiss2_mLc)
        h_D0_ID.Fill(t.Lb_TrueHadron_D0_ID)
        h_D1_ID_2=FillHistoD1(h_D1_ID_2,t.Lb_TrueHadron_D1_ID)
        h_D2_ID.Fill(t.Lb_TrueHadron_D2_ID)
    if t.Lb_TrueHadron_D0_ID==0 and t.Lb_TrueHadron_D1_ID==0 and t.Lb_TrueHadron_D2_ID==0:
        h_E_0.Fill(t.FitVar_El_mLc)
        h_q2_0.Fill(t.FitVar_q2_mLc)
        h_Mmiss2_0.Fill(t.FitVar_Mmiss2_mLc)
        #print(t.Lb_BKGCAT)

leg = TLegend(.60,.60,.85,.80)
leg.SetBorderSize(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(h_E_2b,"2 body decay","L")
leg.AddEntry(h_E_3b,"3 body decay","L")
c = TCanvas('c','c',1500,500)
c.Divide(3,1)
c.cd(1)
h_E_3b.SetLineColor(r.kRed)
h_E_3b.Draw()
h_E_2b.Draw('same')
leg.Draw()
c.cd(2)
h_q2_3b.SetLineColor(r.kRed)
h_q2_3b.Draw()
h_q2_2b.Draw('same')
c.cd(3)
h_Mmiss2_3b.SetLineColor(r.kRed)
h_Mmiss2_3b.Draw()
h_Mmiss2_2b.Draw('same')

c1 = TCanvas('c1','c1',1500,500)
c1.Divide(3,1)
c1.cd(1)
h_E_0.Draw()
c1.cd(2)
h_q2_0.Draw()
c1.cd(3)
h_Mmiss2_0.Draw()

c2 = TCanvas('c2','c2',2100,700)
c2.Divide(3,1)
c2.cd(1)
h_D0_ID.Draw()
c2.cd(2)
h_D1_ID_1.Draw()
h_D1_ID_2.SetLineColor(r.kRed)
h_D1_ID_2.Draw('same')
c2.cd(3)
h_D2_ID.Draw()


print(t.GetEntries(), twobody, threebody)
