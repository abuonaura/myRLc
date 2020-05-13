'''
Author: A. Buonaura
Date: 13th November 2018
Scope: We have observed a weird shape of the Lb_M distribution in same sign data sample. It could be due to the presence of prompt Lambda_c. To verify this hypothesis we plot the Lc_IPCHI2 in the same sign data sample and in the MC sample.


'''

import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas
from ROOT import TMath as m

datadir = '$FILEDIR/'

fd = TFile(datadir+'ControlSamples/Lb_Data_MagUp_sw_noMISID_Kenr.root','READ')
td = fd.Get('DecayTree')

fdf = TFile(datadir+'ControlSamples/Lb_FakeMu_MagUp_reduced_preselected_Kenr_sw.root','READ')
tdf = fdf.Get('DecayTree')

fdss = TFile(datadir+'ControlSamples/Lb_DataSS_MagUp_sw_noMISID_Kenr.root','READ')
tdss = fdss.Get('DecayTree')

fdfss = TFile(datadir+'ControlSamples/Lb_FakeMuSS_MagUp_reduced_preselected_Kenr_sw.root','READ')
tdfss = fdfss.Get('DecayTree')

fmc = TFile(datadir+'MC/Lb_Lcmunu_PID_reduced_preselected.root')
tmc = fmc.Get('DecayTree')

h_IPCHI2d = TH1F('h_IPCHI2d','IP CHI2 data', 100,-5,5) 
h_IPCHI2dss = TH1F('h_IPCHI2dss','IP CHI2 dataSS', 100,-5,5) 
h_IPCHI2df = TH1F('h_IPCHI2df','IP CHI2 data FakeMu', 100,-5,5) 
h_IPCHI2dfss = TH1F('h_IPCHI2dfss','IP CHI2 dataSS FakeMu', 100,-5,5) 
h_IPCHI2mc = TH1F('h_IPCHI2mc','IP CHI2 MC (Lcmunu)', 100,-5,5)


h_LbMMd = TH1F('h_LbMMd','Lb visible mass data', 100,0,10000) 
h_LbMMdss = TH1F('h_LbMMdss','Lb visible mass dataSS', 100,0,10000) 
h_LbMMdf = TH1F('h_LbMMdf','Lb visible mass  data FakeMu', 100,0,10000) 
h_LbMMdfss = TH1F('h_LbMMdfss','Lb visible mass  dataSS FakeMu', 100,0,10000) 
h_LbMMmc = TH1F('h_LbMMmc','Lb visible mass  MC (Lcmunu)', 100,0,10000)

h_LbMmc = TH1F('h_LbMmc','Lb mass  MC (Lcmunu) for Lb_BKGCAT>60', 100,0,10000)
h_LbMdss_cut = TH1F('h_LbMdss_cut','Lb_M evts with IPCHI2>1',100,0,10000)

h_muTrueID = TH1F('h_muTrueID','',8000,-4000,4000)

h_Eld = TH1F('h_Eld','; E_{#mu} (MeV);',20,0,2500)
h_q2d = r.TH1F('h_q2d',';q^{2} (MeV^{2});',20,-2E6,18E6)
h_Mmissd = r.TH1F('h_Mmissd',';Mmiss^{2} (MeV^{2});',20,-2E6,18E6)

h_Eldss = TH1F('h_Eldss','; E_{#mu} (MeV);',20,0,2500)
h_q2dss = r.TH1F('h_q2dss',';q^{2} (MeV^{2});',20,-2E6,18E6)
h_Mmissdss = r.TH1F('h_Mmissdss',';Mmiss^{2} (MeV^{2});',20,-2E6,18E6)

for i in range(td.GetEntries()):
    td.GetEntry(i)
    h_LbMMd.Fill(td.Lb_MM,td.sw_sig)
    h_IPCHI2d.Fill(m.Log10(td.Lc_IPCHI2_OWNPV), td.sw_sig)

for i in range(tdf.GetEntries()):
    tdf.GetEntry(i)
    h_LbMMdf.Fill(tdf.Lb_MM,tdf.sw_sig)
    h_IPCHI2df.Fill(m.Log10(tdf.Lc_IPCHI2_OWNPV), tdf.sw_sig)

for i in range(tdss.GetEntries()):
    tdss.GetEntry(i)
    h_LbMMdss.Fill(tdss.Lb_MM,tdss.sw_sig)
    h_IPCHI2dss.Fill(m.Log10(tdss.Lc_IPCHI2_OWNPV), tdss.sw_sig)

    if (m.Log10(tdss.Lc_IPCHI2_OWNPV)>1):
        h_LbMdss_cut.Fill(tdss.Lb_M, tdss.sw_sig)

for i in range(tdfss.GetEntries()):
    tdfss.GetEntry(i)
    h_LbMMdfss.Fill(tdfss.Lb_MM,tdfss.sw_sig)
    h_IPCHI2dfss.Fill(m.Log10(tdfss.Lc_IPCHI2_OWNPV), tdfss.sw_sig)

for j in range(tmc.GetEntries()):
    tmc.GetEntry(j)
    if(tmc.Lc_BKGCAT<30):
        if (tmc.Lb_ISOLATION_BDT> 0.35 and tmc.Lb_ISOLATION_BDT2>0.2)and((tmc.Lb_ISOLATION_PIDK>4.and(tmc.Lb_ISOLATION_CHARGE==tmc.mu_ID/13 or (tmc.Lb_ISOLATION_CHARGE==-tmc.mu_ID/13 and tmc.Lb_ISOLATION_PIDp - tmc.Lb_ISOLATION_PIDK<0.))) or (tmc.Lb_ISOLATION_PIDK2>4.and(tmc.Lb_ISOLATION_CHARGE2==tmc.mu_ID/13 or (tmc.Lb_ISOLATION_CHARGE2==-tmc.mu_ID/13 and tmc.Lb_ISOLATION_PIDp2 - tmc.Lb_ISOLATION_PIDK2<0.)))):
            h_LbMMmc.Fill(tmc.Lb_MM)
            h_IPCHI2mc.Fill(m.Log10(tmc.Lc_IPCHI2_OWNPV))
            if(tmc.Lb_BKGCAT>=60):
                h_LbMmc.Fill(tmc.Lb_M)
                h_muTrueID.Fill(tmc.mu_TRUEID)

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

h_IPCHI2d=ScaleHisto(h_IPCHI2d,1)
h_IPCHI2df=ScaleHisto(h_IPCHI2df,1)
h_IPCHI2dss=ScaleHisto(h_IPCHI2dss,1)
h_IPCHI2dfss=ScaleHisto(h_IPCHI2dfss,1)
h_IPCHI2mc=ScaleHisto(h_IPCHI2mc,1)

h_LbMMd = ScaleHisto(h_LbMMd,1)
h_LbMMdf = ScaleHisto(h_LbMMdf,1)
h_LbMMdss = ScaleHisto(h_LbMMdss,1)
h_LbMMdfss = ScaleHisto(h_LbMMdfss,1)
h_LbMMmc = ScaleHisto(h_LbMMmc,1)



legend=r.TLegend(0.6,0.6,0.9,0.9)
legend.AddEntry(h_LbMMd,'Data', "l")
legend.AddEntry(h_LbMMdf,'Data FakeMu', "l")
legend.AddEntry(h_LbMMdss,'Data SameSign', "l")
legend.AddEntry(h_LbMMdfss,'Data SameSign FakeMu', "l")
legend.AddEntry(h_LbMMmc,'MC (#Lambda_{b}->#Lambda_{c}#mu#nu)', "l")

c = TCanvas('c','',500,500)
h_IPCHI2dss.SetLineColor(r.kRed)
h_IPCHI2d.SetLineColor(r.kBlack)
h_IPCHI2df.SetLineColor(r.kGreen)
h_IPCHI2dfss.SetLineColor(r.kViolet)
h_IPCHI2mc.Draw('hist')
h_IPCHI2df.Draw('hist sames')
h_IPCHI2dfss.Draw('hist sames')
h_IPCHI2dss.Draw('hist sames')
h_IPCHI2d.Draw('hist sames')
legend.Draw()

c1 = TCanvas('c1','',500,500)
h_LbMMmc.GetYaxis().SetRangeUser(0,0.13)
h_LbMMdss.SetLineColor(r.kRed)
h_LbMMd.SetLineColor(r.kBlack)
h_LbMMdf.SetLineColor(r.kGreen)
h_LbMMdfss.SetLineColor(r.kViolet)
h_LbMMmc.Draw('hist' )
h_LbMMdf.Draw('hist sames')
h_LbMMdfss.Draw('hist sames')
h_LbMMdss.Draw('hist sames')
h_LbMMd.Draw('hist sames')
legend.Draw()

c2 = TCanvas('c2','',500,500)
h_LbMmc.Draw('hist')

c3 = TCanvas('c3','',500,500)
h_LbMdss_cut.Draw('hist')

c4 = TCanvas('c4','',500,500)
h_muTrueID.Draw('hist')

