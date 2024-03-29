import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas
import numpy as np

ISOBDTcut =0.35
ISOBDT2cut=0.2
#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)

def CheckIfIsKenriched(t,evt):
    t.GetEntry(evt)
    BDT = t.Lb_ISOLATION_BDT
    BDT2 = t.Lb_ISOLATION_BDT2
    PIDp = t.Lb_ISOLATION_PIDp
    PIDK = t.Lb_ISOLATION_PIDK
    Ch = t.Lb_ISOLATION_CHARGE
    Type = t.Lb_ISOLATION_Type
    NNghost = t.Lb_ISOLATION_NNghost
    PIDp2 = t.Lb_ISOLATION_PIDp2
    PIDK2 =t.Lb_ISOLATION_PIDK2
    Ch2 = t.Lb_ISOLATION_CHARGE2
    Type2 = t.Lb_ISOLATION_Type2
    NNghost2 = t.Lb_ISOLATION_NNghost2
    muCh = -t.mu_ID/13
    isK1=0
    isK2=0
    #print(BDT,BDT2,PIDp,PIDK,Ch,Type,NNghost,PIDp2,PIDK2,Type2,Ch2,NNghost2,muCh)

    if BDT>ISOBDTcut and BDT2>ISOBDT2cut:
        m1 = m_pi
        m2 = m_pi
	#I retrieve the momenta of the 2 anti-isolated particles:
        p1 = np.matrix([t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ])
        p2 = np.matrix([t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2])
        #------------------------
        #I will do in the following several mass hypothesis for these 2 particles to make
        #the Kenriched sample as cleaner as possible
        #------------------------
        #I want to discard events coming from LcStar ->Lcpipi
        #-> I assume that the 2 anti-isolated particles are 2 pions
        E1pi = r.TMath.Sqrt(p1*p1.transpose() + m_pi*m_pi)
        E2pi = r.TMath.Sqrt(p2*p2.transpose() + m_pi*m_pi)
        #-> I retrieve the momentum of the Lambda_c
        pLc = np.matrix([t.Lc_PX, t.Lc_PY, t.Lc_PZ])
        ELc = r.TMath.Sqrt(pLc*pLc.transpose() + m_Lc*m_Lc)
        #-> I compute the tot. energy of the 3 particle system
        ELc12 = E1pi+E2pi+ELc
        pLc12 = p1+p2+pLc
        mLc12 = r.TMath.Sqrt(ELc12*ELc12 - pLc12*pLc12.transpose())
        #    h_Lcpipi.Fill(mLc12)
        #------------------------
        #I want to discard events with a total mass above Lbmass
        #This time I check which of the 2 (or both) looks more like a K
        #I consider the other particle to be a pion
        if PIDK>4.:
            if (Ch==-muCh or (Ch==muCh and (PIDp-PIDK)<0.)):
                m1 = m_K
                isK1=1
            if (Ch==muCh and PIDp-PIDK>0.):#I assume it is a pion
                m1 = m_pi
        if PIDK2>4.:
            if (Ch2==-muCh or (Ch2==muCh and (PIDp2 - PIDK2)<0.)):
                #I also assume this particle is a K
                m2 = m_K
                isK2=1
            if (Ch2==muCh and PIDp2-PIDK2>0.):
                m2 = m_pi
        #-> I compute the energy of the 2 particles:
        E1 = r.TMath.Sqrt(p1*p1.transpose() + m1*m1)
        E2 = r.TMath.Sqrt(p2*p2.transpose() + m2*m2)
        #-> I retrieve the momentum of the muon
        pmu = np.matrix([t.mu_PX, t.mu_PY, t.mu_PZ])
        Emu = r.TMath.Sqrt(pmu*pmu.transpose() + m_mu*m_mu)
        #-> I compute the total momentum and energy of the 4 particles:
        pTOT = p1 + p2 + pmu + pLc
        ETOT = E1 + E2 + Emu + ELc
        #-> I evaluate the invariant mass
        mTOT = r.TMath.Sqrt(ETOT*ETOT - pTOT*pTOT.transpose())

        #print(m2,m1,isK1, isK2,mTOT,mLc12)
	#I save only those events with mLc12>2700 and with mTOT<5620
        if mLc12>2700 and mTOT<5620 and (isK1==1 or isK2==1):
            isKenriched=1
        else:
            isKenriched=0
    else:
        isKenriched=0
    
    return isKenriched
		
def CheckIfIsIsolated(t,evt):
    t.GetEntry(evt)

    BDT = t.Lb_ISOLATION_BDT
    if BDT<ISOBDTcut:
        isIsolated=1
    else:
        isIsolated=0
    return isIsolated


