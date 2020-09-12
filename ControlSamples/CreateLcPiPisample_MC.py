import ROOT as r

import os,sys,getopt,time
import numpy as np
import argparse
sys.path.append('../Preselection/')
from AddSweights import *



datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
polarities = ['MagUp','MagDown']

ISOBDTcut = 0.35
ISOBDT2cut = 0.2

#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)

def SearchDecay(t):
    h_mLcSt = r.TH1F('h_mLcSt',';m(GeV/c^{2});',200,2200,5000)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if i%100000==0:
            print(i)
        if t.FinalSel==True:
            weight = t.Event_PIDCalibEffWeight
            BDT = t.Lb_ISOLATION_BDT
            PIDp = t.Lb_ISOLATION_PIDp
            PIDK = t.Lb_ISOLATION_PIDK
            Ch   = t.Lb_ISOLATION_CHARGE
            Type = t.Lb_ISOLATION_Type
            NNghost = t.Lb_ISOLATION_NNghost
            BDT2 = t.Lb_ISOLATION_BDT2
            PIDp2 = t.Lb_ISOLATION_PIDp2
            PIDK2 = t.Lb_ISOLATION_PIDK2
            Ch2 = t.Lb_ISOLATION_CHARGE2
            Type2 = t.Lb_ISOLATION_Type2
            NNghost2 = t.Lb_ISOLATION_NNghost2
            BDT3 = t.Lb_ISOLATION_BDT3
            PIDp3 = t.Lb_ISOLATION_PIDp3
            PIDK3 = t.Lb_ISOLATION_PIDK3
            Ch3   = t.Lb_ISOLATION_CHARGE3
            Type3 = t.Lb_ISOLATION_Type3
            NNghost3 = t.Lb_ISOLATION_NNghost3
            if BDT>0.35 and BDT2>0.2:
                #I require that the 2 anti-isolated particles have opposite charge
                if Ch==-Ch2:
                    #I assign them the mass of the pion
                    m1 = m_pi
                    m2 = m_pi
                    #I retrieve the momenta of the 2 anti-isolated particles:
                    p1 = np.matrix([t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ])
                    p2 = np.matrix([t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2])
                    #I compute the energy of the 2 particles:
                    E1 = r.TMath.Sqrt(p1*p1.transpose() + m1*m1)
                    E2 = r.TMath.Sqrt(p2*p2.transpose() + m2*m2)
                    #I retrieve the momentum of the Lambda_c
                    pLc = np.matrix([t.Lc_PX, t.Lc_PY, t.Lc_PZ])
                    ELc = r.TMath.Sqrt(pLc*pLc.transpose() + m_Lc*m_Lc)

                    #I compute the total momentum and energy of this 3 particles system
                    pLcSt = p1+p2+pLc
                    ELcSt = E1+E2+ELc
                    #and compute the invariant mass:
                    mLcSt = r.TMath.Sqrt(ELcSt*ELcSt - pLcSt*pLcSt.transpose())
                    h_mLcSt.Fill(mLcSt, weight)
    h_mLcSt.SetDirectory(0)
    return h_mLcSt

fname = datadir+'Lb_Lcmunu_MagUp_full.root'
preselfname = fname[0:-5]+'_preselectionVars.root'
pidcalibFile = fname[0:-5]+'_PIDCalib.root'

f = r.TFile(fname,'READ')
t = f.Get('tupleout/DecayTree')
fpresel = r.TFile(preselfname,'READ')
tpresel = fpresel.Get('DecayTree')
t.AddFriend(tpresel)
fpidcalib = r.TFile(pidcalibFile)
tcalib = fpidcalib.Get('tupleout/DecayTree')
t.AddFriend(tcalib)

h_mLcSt = SearchDecay(t)

'''
#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)
#Histograms
h_E1pi = r.TH1F('h_E1pi','',100,0,50000)
h_E2pi = r.TH1F('h_E2pi','',100,0,50000)
h_ELc = r.TH1F('h_ELc','',100,0,50000)
h_mLc12 = r.TH1F('h_mLc12','',100,2000,5000)
h_E1pi2 = r.TH1F('h_E1pi2','',100,0,50000)
h_E2pi2= r.TH1F('h_E2pi2','',100,0,50000)
h_ELc2 = r.TH1F('h_ELc2','',100,0,50000)
h_mLc12_2 = r.TH1F('h_mLc12_2','',100,2000,5000)
for i in range(100000):
    t.GetEntry(i)
    if i%100000==0:
        print(i)
    #Method 1
    E1pi = r.TMath.Sqrt(r.TMath.Power(t.Lb_ISOLATION_PX,2)+r.TMath.Power(t.Lb_ISOLATION_PY,2)+r.TMath.Power(t.Lb_ISOLATION_PZ,2)+r.TMath.Power(m_pi,2))
    E2pi = r.TMath.Sqrt(r.TMath.Power(t.Lb_ISOLATION_PX2,2)+r.TMath.Power(t.Lb_ISOLATION_PY2,2)+r.TMath.Power(t.Lb_ISOLATION_PZ2,2)+r.TMath.Power(m_pi,2))
    ELc = r.TMath.Sqrt(r.TMath.Power(t.Lc_PX,2)+r.TMath.Power(t.Lc_PY,2)+r.TMath.Power(t.Lc_PZ,2)+m_Lc*m_Lc)
    pLc12_x = t.Lc_PX + t.Lb_ISOLATION_PX + t.Lb_ISOLATION_PX2
    pLc12_y = t.Lc_PY+t.Lb_ISOLATION_PY+t.Lb_ISOLATION_PY2
    pLc12_z = t.Lc_PZ+t.Lb_ISOLATION_PZ+t.Lb_ISOLATION_PZ2
    mLc12 = r.TMath.Sqrt(r.TMath.Power(ELc+E1pi+E2pi,2)-(r.TMath.Power(pLc12_x,2)+r.TMath.Power(pLc12_y,2)+r.TMath.Power(pLc12_z,2)))
    h_E1pi.Fill(E1pi)
    h_E2pi.Fill(E2pi)
    h_ELc.Fill(ELc)
    h_mLc12.Fill(mLc12)
    #Method2
    p1 = np.matrix([t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ])
    p2 = np.matrix([t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2])
    E1pi2 = r.TMath.Sqrt(p1*p1.transpose() + m_pi*m_pi)
    E2pi2 = r.TMath.Sqrt(p2*p2.transpose() + m_pi*m_pi)
    pLc = np.matrix([t.Lc_PX, t.Lc_PY, t.Lc_PZ])
    ELc2 = r.TMath.Sqrt(pLc*pLc.transpose() + m_Lc*m_Lc)
    h_E1pi2.Fill(E1pi2)
    h_E2pi2.Fill(E2pi2)
    h_ELc2.Fill(ELc2)
    pLcSt = p1+p2+pLc
    ELcSt = E1pi2+E2pi2+ELc
    #and compute the invariant mass:
    mLc12_2 = r.TMath.Sqrt(ELcSt*ELcSt - pLcSt*pLcSt.transpose())
    h_mLc12_2.Fill(mLc12_2)

c = r.TCanvas('c','c',1000,500)
c.Divide(2,1)
c.cd(1)
h_E1pi.Draw()
h_E1pi2.SetLineColor(r.kRed)
h_E1pi2.Draw('sames')
c.cd(2)
h_E2pi.Draw()
h_E2pi2.SetLineColor(r.kRed)
h_E2pi2.Draw('sames')

c1 = r.TCanvas('c1','c1')
h_ELc.Draw()
h_ELc2.SetLineColor(r.kRed)
h_ELc2.Draw('sames')

c2 = r.TCanvas('c2','c2')
h_mLc12.Draw()
h_mLc12_2.SetLineColor(r.kRed)
h_mLc12_2.Draw('sames')
'''
