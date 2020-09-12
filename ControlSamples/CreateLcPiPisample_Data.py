import ROOT as r

import os,sys,getopt,time
import numpy as np
import argparse
sys.path.append('../Preselection/')
from AddSweights import *


parser = argparse.ArgumentParser(description='Split samples accoriding to isolation requests + add sweights')
parser.add_argument('datatype',choices=['Data','FakeMu','FakeMuSS','DataSS','all'], help = 'which data sample we want to run on')
parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on')

args = parser.parse_args()

datatype =args.datatype
polarity=args.polarity

datadir = '/disk/lhcb_data2/RLcMuonic2016/Data/'
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
        if i%1000000==0:
            print(i)
        if t.FinalSel==True and t.isKenriched==False:
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
                    h_mLcSt.Fill(mLcSt)
    h_mLcSt.SetDirectory(0)
    return h_mLcSt

def PlotSweightedLcpipi(fname):
    f = r.TFile(fname, 'READ')
    t = f.Get('DecayTree')
    h_mLcSt = r.TH1F('h_mLcSt',';m(GeV/c^{2});',200,2200,5000)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if i%100000==0:
            print(i)
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
        h_mLcSt.Fill(mLcSt,t.sw_sig)
    h_mLcSt.SetDirectory(0)
    return h_mLcSt




def SplitSample(fname, t):
    ofname = fname[0:-5]+'_Lcpipi.root'
    of = r.TFile(ofname,'recreate')
    ot = r.TTree('DecayTree','DecayTree')
    print(' >>> Copying tree')
    ot = t.CloneTree(0)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if i%1000000==0:
            print(i)
        if t.FinalSel==True and t.isKenriched==False:
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
                    if mLcSt<2700:
                        ot.Fill()
    of.Write()
    of.Close()
    return ofname

if __name__== "__main__":
    if datatype=='all':
        datatypes = ['Data','FakeMu','FakeMuSS','DataSS']
    else:
        datatypes=[datatype]
    if polarity!='all':
        polarities=[polarity]

    for dtype in datatypes:
        print('>>>>   Processing %s sample' %(dtype))
        for pol in polarities:
            fname = datadir+'Lb_'+dtype+'_'+pol+'.root'
            preselfname = fname[0:-5]+'_Df_preselection.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            fpresel = r.TFile(preselfname,'READ')
            tpresel = fpresel.Get('DecayTree')
            t.AddFriend(tpresel)

            h_mLcSt = SearchDecay(t)
            h_mLcSt.Draw()
            ''' 
            ofname = SplitSample(fname, t)
            #ofname = fname[0:-5]+'_Lcpipi.root'
            swfname =ComputeSweights(ofname, dtype, pol,'Lcpipi')
            tname = 'DecayTree'
            variable = 'Lc_M'
            splotVariable(swfname, tname, dtype, pol,'Lcpipi',variable, 50, 2230, 2330, 'm_{#Lambda_{c}}')
            swfname = fname[0:-5]+'_Lcpipi_sw.root'
            hmLcSt = PlotSweightedLcpipi(swfname)
            '''




