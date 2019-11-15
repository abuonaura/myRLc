import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas, TMath, TLegend
import numpy as np

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'

fname = datadir+'MC/Lb_LcDs_PID_reduced_preselected.root'

f = TFile(fname,'READ')
t = f.Get('DecayTree')

#Masses pi, K, p
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)


h_1mu = r.TH1F('h_1mu',';m (MeV);',100,200,2000)
h_2mu = r.TH1F('h_2mu',';m (MeV);',100,200,2000)
h_12 = r.TH1F('h_12',';m (MeV);',100,100,3000)
h_12mu = r.TH1F('h_12mu',';m (MeV);',200,200,5000)
h_TOT = r.TH1F('h_TOT',';m (MeV);',200,3000,8000)
h_LbM = r.TH1F('h_LbM',';m (MeV);',200,1000,7000)

for i in range(t.GetEntries()):
    t.GetEntry(i)
    PIDp = t.Lb_ISOLATION_PIDp
    PIDK = t.Lb_ISOLATION_PIDK
    Ch = t.Lb_ISOLATION_CHARGE
    Type = t.Lb_ISOLATION_Type
    NNghost = t.Lb_ISOLATION_NNghost
    PIDp2 = t.Lb_ISOLATION_PIDp2
    PIDK2 = t.Lb_ISOLATION_PIDK2
    Ch2 = t.Lb_ISOLATION_CHARGE2
    Type2 = t.Lb_ISOLATION_Type2
    NNghost2 = t.Lb_ISOLATION_NNghost2
    muCh = -t.mu_ID/13
    #truth matching
    if Type==3 and NNghost<0.2 and Type2==3 and NNghost2<0.2 and t.Lb_BKGCAT<50. and t.Lc_BKGCAT<30.:
        if t.Lb_ISOLATION_BDT>0.35 and t.Lb_ISOLATION_BDT2>0.2:
            #I retrieve the momenta of the 2 anti-isolated particles:
            p1 = np.matrix([t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ])
            p2 = np.matrix([t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2])
            #I see which mass hypothesis is more suitable for the 2 anti isolated particles
            if PIDK>4.:
                if (Ch==-muCh or (Ch==muCh and (PIDp-PIDK)<0.)):
                    #I also assume this particle is a K
                    m1 = m_K 
                if (Ch==muCh and PIDp-PIDK>0.):
                    #In this case I assume it is a proton
                    m1 = m_p
            else: #if it does not satisfy the hypothesis for being a K, I assume it is a pi
                m1 = m_pi
            #and I check which masses I can assume for the second
            if PIDK2>4.: 
                if (Ch2==-muCh or (Ch2==muCh and (PIDp2 - PIDK2)<0.)):
                    #I also assume this particle is a K
                    m2 = m_K
                if (Ch2==muCh and PIDp2-PIDK2>0.):
                    #In this case I assume it is a proton
                    m2 = m_p
            else: #if it does not satisfy the hypothesis for being a K, I assume it is a pi
                m2 = m_pi

            #I compute the energy of the 2 particles:
            E1 = r.TMath.Sqrt(p1*p1.transpose() + m1*m1)
            E2 = r.TMath.Sqrt(p2*p2.transpose() + m2*m2)

            #I retrieve the momentum of the muon
            pmu = np.matrix([t.mu_PX, t.mu_PY, t.mu_PZ])
            Emu = r.TMath.Sqrt(pmu*pmu.transpose() + m_mu*m_mu)

            #I evaulate the total energy and momentum of the following combinations: 1+mu, 2+mu, 1+2, 1+2+mu
            E1mu = E1+Emu
            p1mu = p1+pmu

            E2mu = E2+Emu
            p2mu = p2+pmu

            E12 = E1+E2
            p12 = p1+p1

            E12mu = E1+E2+Emu
            p12mu = p1+p2+pmu

            #and I retrieve the invariant mass for the 4 hypothesis
            m1mu = r.TMath.Sqrt(E1mu*E1mu - p1mu*p1mu.transpose())
            m2mu = r.TMath.Sqrt(E2mu*E2mu - p2mu*p2mu.transpose())
            m12 = r.TMath.Sqrt(E12*E12 - p12*p12.transpose())
            m12mu = r.TMath.Sqrt(E12mu*E12mu - p12mu*p12mu.transpose())

            h_1mu.Fill(m1mu)
            h_2mu.Fill(m2mu)
            h_12.Fill(m12)
            h_12mu.Fill(m12mu)

            #I retrieve the momentum of the Lambda_c
            pLc = np.matrix([t.Lc_PX, t.Lc_PY, t.Lc_PZ])
            ELc = r.TMath.Sqrt(pLc*pLc.transpose() + m_Lc*m_Lc)

            #I compute the total momentum and energy of the 4 particles:
            pTOT = p12mu + pLc
            ETOT = E12mu + ELc
            #I evaluate the invariant mass
            mTOT = r.TMath.Sqrt(ETOT*ETOT - pTOT*pTOT.transpose())
            h_TOT.Fill(mTOT)

            h_LbM.Fill(t.Lb_M)

c = TCanvas('c','1+mu',500,500)
h_1mu.Draw()

c1 = TCanvas('c1','2+mu',500,500)
h_2mu.Draw()

c3 = TCanvas('c3','1+2',500,500)
h_12.Draw()

c4 = TCanvas('c4','1+2+mu',500,500)
h_12mu.Draw()                    

c5 = TCanvas('c5','tot',500,500)
h_TOT.Draw()
            
c6 = TCanvas('c6','LbM',500,500)
h_LbM.Draw()
            
            

            
            
            
            #if ((t.Lb_ISOLATION_PIDK>4. and (t.Lb_ISOLATION_CHARGE=-t.mu_ID/13.) or (t.Lb_ISOLATION_CHARGE=t.mu_ID/13. and t.Lb_ISOLATION_PIDp-t.Lb_ISOLATION_PIDK<0.)) or (t.Lb_ISOLATION_PIDK2>4. and ((t.Lb_ISOLATION_CHARGE2=-t.mu_ID/13.) or (t.Lb_ISOLATION_CHARGE2=t.mu_ID/13. and t.Lb_ISOLATION_PIDp2-t.Lb_ISOLATION_PIDK2<0.)))):


