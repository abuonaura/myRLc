import ROOT as r
from ROOT import TMath as mt
import numpy as np

#Starting variables
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)
ISOBDTcut=0.35
ISOBDT2cut=0.2
Kcodes = [130,310,311,321,313,323,10311,10321,10313,10323,20313,20323,30313,30323,10315,10325,20315,20325,317,327,319,329]

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h


#Upload files to read
datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
samples = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
samples = ['Lcmunu','LcDs']
#samples = ['LcDs']

for sample in samples:
    print('--------------- '+sample+' ---------------')
    f = r.TFile(datadir+'Lb_'+sample+'_MagUp.root','READ')
    f1 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_preselectionVars.root')
    f2 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_PIDCalib.root')

    #Get Tree + tree with preselection variables
    t = f.Get('tupleout/DecayTree')
    t1 = f1.Get('DecayTree')
    t2 = f2.Get('tupleout/DecayTree')

    t.AddFriend(t1)
    t.AddFriend(t2)
    #NOTE: The number of entries passing the selection is not the number of events, because it has to be reweighted for PID calib
    print('Number of unpreselected entries: ', t.GetEntries())
    print('Number of enties passing the final selection: ',t.Draw('Lc_M','FinalSel==1','goff'))

    #Switch off all branches to speed up
    t.SetBranchStatus('*',0)
    #List of branches to turn on
    blist = ['FinalSel','isKenriched','Lb_ISOLATION_TruePID','Lb_ISOLATION_TrueMotherPID','Lb_ISOLATION_TruePID2','Lb_ISOLATION_TrueMotherPID2','Lb_ISOLATION_BDT','Lb_ISOLATION_BDT2','Lb_ISOLATION_BDT3','Lb_ISOLATION_PIDK','Lb_ISOLATION_PIDp','Lb_ISOLATION_PIDK2','Lb_ISOLATION_PIDp2','Lb_ISOLATION_PIDK3','Lb_ISOLATION_PIDp3','Lb_ISOLATION_Type', 'Lb_ISOLATION_Type2','Lb_ISOLATION_Type3','Lb_ISOLATION_NNghost', 'Lb_ISOLATION_NNghost2','Lb_ISOLATION_NNghost3', 'Lb_ISOLATION_CHARGE', 'Lb_ISOLATION_CHARGE2','Lb_ISOLATION_CHARGE3','Event_PIDCalibEffWeight','mu_ID','Lb_ISOLATION_PX','Lb_ISOLATION_PX2','Lb_ISOLATION_PY','Lb_ISOLATION_PY2','Lb_ISOLATION_PZ','Lb_ISOLATION_PZ2','Lc_PX','Lc_PY','Lc_PZ','mu_PX','mu_PY','mu_PZ','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc']
    for b in blist:
        t.SetBranchStatus(b,1)

    #histos to be filled with the templates
    h_std = r.TH3F('h_std','Standard Kenriched Selection '+sample+';q^{2};E_{#mu};M_{miss}^{2}',4,-2,14,10,0,2600,10,-2,14)
    h_new = r.TH3F('h_new','New Kenriched Selection '+sample+';q^{2};E_{#mu};M_{miss}^{2}',4,-2,14,10,0,2600,10,-2,14)

    h_TrueID = r.TH1F('h_TrueID','TrueID after asking PIDK>4; |TrueID|',2500,0,2500)
    h_TrueID_1 = r.TH1F('h_TrueID_1','TrueID after asking PIDK>4 and additional K selections; |TrueID|',2500,0,2500)
    h_TrueID2 = r.TH1F('h_TrueID2','TrueID2 after asking PIDK2>4; |TrueID2|',2500,0,2500)
    h_TrueID2_1 = r.TH1F('h_TrueID2_1','TrueID2 after asking PIDK2>4 and additional K selections; |TrueID2|',2500,0,2500)

    h_PIDK = r.TH1F('h_PIDK','PIDK for K particles',50,-100,100)
    h_PIDK_1 = r.TH1F('h_PIDK_1','PIDK for K particles',50,-100,100)
    h_PIDK2 = r.TH1F('h_PIDK2','PIDK2 for K particles',50,-100,100)
    h_PIDK2_1 = r.TH1F('h_PIDK2_1','PIDK2 for K particles',50,-100,100)
    h_pK = r.TH3F('h_pK','Momentum of True K particles; p_{x}; p_{y}; p_{z}',50,-3E3,3E3,50,-3E3,3E3,50,0,2E5)
    h_pK_1 = r.TH3F('h_pK_1','Momentum of particles with PIDK>4; p_{x}; p_{y}; p_{z}',50,-3E3,3E3,50,-3E3,3E3,50,0,2E5)
    h_pK2 = r.TH3F('h_pK2','Momentum of True K particles (2); p_{x}; p_{y}; p_{z}',50,-3E3,3E3,50,-3E3,3E3,50,0,2E5)
    h_pK2_1 = r.TH3F('h_pK2_1','Momentum of particles with PIDK2>4; p_{x}; p_{y}; p_{z}',50,-3E3,3E3,50,-3E3,3E3,50,0,2E5)

    h_PIDKp = r.TH2F('h_PIDKp','PIDK Vs p; PIDK; p',50,-10,10,50,0,5.E4)
    h_PIDK2p2= r.TH2F('h_PIDK2p2','PIDK2 Vs p2; PIDK2; p2',50,-10,10,50,0,5E4)
    
    h_PIDdiff = r.TH1F('h_PIDdiff','PIDp-PIDK for true K',100,-500,500)
    h_PIDdiff_1 = r.TH1F('h_PIDdiff_1','PIDp - PIDK for true K with PIDK>4',100,-500,500)
    h_PIDdiff2 = r.TH1F('h_PIDdiff2','PIDp-PIDK for true K (2)',100,-500,500)
    h_PIDdiff2_1 = r.TH1F('h_PIDdiff2_1','PIDp - PIDK for true K (2) with PIDK>4',100,-500,500)

    h_mLc12 = r.TH1F('h_mLc12','',100,2000,5000)
    h_mTOT_old = r.TH1F('h_mLcTOT_old','',100,2000,7000)
    h_mTOT_new = r.TH1F('h_mLcTOT_new','',100,2000,7000)
    h_E1pi = r.TH1F('h_E1pi','',100,0,10000)
    h_E2pi = r.TH1F('h_E2pi','',100,0,10000)
    h_ELc = r.TH1F('h_ELc','',100,0,1E6)
    h_m1 = r.TH1F('h_m1','',100,0,1000)
    h_m2 = r.TH1F('h_m2','',100,0,1000)
    h_muCh = r.TH1F('h_muCh','',3,-1,1)


    #Loop over entries
    for i in range(t.GetEntries()):
    #for i in range(500000):
        t.GetEntry(i)
        #discard events not passing all the selection requirements
        if t.FinalSel!=1:
            continue
        else:
            if t.isKenriched==1:
                h_std.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,t.Event_PIDCalibEffWeight)

            #Variables necessary for studying new Kenriched selection
            isKenriched_new=0
            muCh = -t.mu_ID/13
            h_muCh.Fill(muCh)
            BDT = t.Lb_ISOLATION_BDT
            PIDp = t.Lb_ISOLATION_PIDp
            PIDK = t.Lb_ISOLATION_PIDK
            Ch   = t.Lb_ISOLATION_CHARGE
            Type = t.Lb_ISOLATION_Type
            TrueID = t.Lb_ISOLATION_TruePID
            TrueMID = t.Lb_ISOLATION_TrueMotherPID
            NNghost = t.Lb_ISOLATION_NNghost
            BDT2 = t.Lb_ISOLATION_BDT2
            PIDp2 = t.Lb_ISOLATION_PIDp2
            PIDK2 = t.Lb_ISOLATION_PIDK2
            Ch2 = t.Lb_ISOLATION_CHARGE2
            Type2 = t.Lb_ISOLATION_Type2
            TrueID2 = t.Lb_ISOLATION_TruePID2
            TrueMID2 = t.Lb_ISOLATION_TrueMotherPID2
            NNghost2 = t.Lb_ISOLATION_NNghost2
            isK1=0
            isK2=0
            E1pi,E2pi,ELc,mLc12=0.,0.,0.,0.
            if Type==3 and NNghost<0.2 and Type2==3 and NNghost2<0.2:
                if BDT>ISOBDTcut and BDT2>ISOBDT2cut:
            
                    #---------- KEnriched selection --------------
                    #I retrieve the momenta of the 2 anti-isolated particles:
                    p1 = np.matrix([t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ])
                    p2 = np.matrix([t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2])
                    h_PIDKp.Fill(PIDK,mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX,2)+mt.Power(t.Lb_ISOLATION_PY,2)+mt.Power(t.Lb_ISOLATION_PZ,2)))
                    h_PIDK2p2.Fill(PIDK2,mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX2,2)+mt.Power(t.Lb_ISOLATION_PY2,2)+mt.Power(t.Lb_ISOLATION_PZ2,2)))
                    #------------------------
                    #I will do in the following several mass hypothesis for these 2 particles to make
                    #the Kenriched sample as cleaner as possible
                    #------------------------
                    #I want to discard events coming from LcStar ->Lcpipi
                    #-> I assume that the 2 anti-isolated particles are 2 pions
                    E1pi = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX,2)+mt.Power(t.Lb_ISOLATION_PY,2)+mt.Power(t.Lb_ISOLATION_PZ,2)+m_pi*m_pi)
                    E2pi = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX2,2)+mt.Power(t.Lb_ISOLATION_PY2,2)+mt.Power(t.Lb_ISOLATION_PZ2,2)+m_pi*m_pi)
                    ELc = mt.Sqrt(mt.Power(t.Lc_PX,2)+mt.Power(t.Lc_PY,2)+mt.Power(t.Lc_PZ,2)+m_Lc*m_Lc)
                    #-> I compute the tot. energy of the 3 particle system
                    pLc12_x = t.Lc_PX+t.Lb_ISOLATION_PX+t.Lb_ISOLATION_PX2
                    pLc12_y = t.Lc_PY+t.Lb_ISOLATION_PY+t.Lb_ISOLATION_PY2
                    pLc12_z = t.Lc_PZ+t.Lb_ISOLATION_PZ+t.Lb_ISOLATION_PZ2
                    mLc12 = mt.Sqrt(mt.Power(ELc+E1pi+E2pi,2)-(mt.Power(pLc12_x,2)+mt.Power(pLc12_y,2)+mt.Power(pLc12_z,2)))
                    #I check if one of the two particles is a K otherwise otherwise I assign the pi mass:
                    m1 = m_pi
                    m1old=m_pi
                    
                    
                    if PIDK>4.:
                        h_TrueID.Fill(r.TMath.Abs(TrueID))
                        if ((Ch==-muCh) or (Ch==muCh and (PIDp-PIDK)<=0.)):
                            m1old=m_K
                            #isK1=1
                            h_TrueID_1.Fill(r.TMath.Abs(TrueID))
                            h_pK_1.Fill(t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ)
                        if (Ch==muCh and (PIDp-PIDK)>0.):
                            m1old=m_p
                        
                    m2 = m_pi
                    m2old = m_pi
                    if PIDK2>4:
                        h_TrueID2.Fill(r.TMath.Abs(TrueID2))
                        if (Ch2==-muCh or (Ch2==muCh and (PIDp2-PIDK2)<=0.)):
                            m2old = m_K
                            #isK2=1
                            h_TrueID2_1.Fill(r.TMath.Abs(TrueID2))
                            h_pK2_1.Fill(t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2)
                        if (Ch2==muCh and (PIDp2-PIDK2)>0.):
                            m2old=m_p

                    if r.TMath.Abs(TrueID) in Kcodes:
                        m1 = m_K
                        isK1=1
                        h_PIDK.Fill(PIDK)
                        h_pK.Fill(t.Lb_ISOLATION_PX, t.Lb_ISOLATION_PY, t.Lb_ISOLATION_PZ)
                        h_PIDdiff.Fill(PIDp-PIDK)
                        if PIDK>4.:
                            h_PIDK_1.Fill(PIDK)
                            h_PIDdiff_1.Fill(PIDp-PIDK)
                    else:
                        if r.TMath.Abs(TrueID)==2212:
                            m1 = m_p
                        else: 
                            m1 = m_pi
                    if r.TMath.Abs(TrueID2) in Kcodes:
                        m2 = m_K
                        h_PIDK2.Fill(PIDK2)
                        h_pK2.Fill(t.Lb_ISOLATION_PX2, t.Lb_ISOLATION_PY2, t.Lb_ISOLATION_PZ2)
                        h_PIDdiff2.Fill(PIDp2-PIDK2)
                        if PIDK2>4.:
                            h_PIDK2_1.Fill(PIDK2)
                            h_PIDdiff2_1.Fill(PIDp2-PIDK2)
                    else:
                        if r.TMath.Abs(TrueID2)==2212:
                            m2 = m_p
                        else: 
                            m2= m_pi

                    E1 = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX,2)+mt.Power(t.Lb_ISOLATION_PY,2)+mt.Power(t.Lb_ISOLATION_PZ,2)+mt.Power(m1,2))
                    E2 = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX2,2)+mt.Power(t.Lb_ISOLATION_PY2,2)+mt.Power(t.Lb_ISOLATION_PZ2,2)+mt.Power(m2,2))
                    
                    E1old = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX,2)+mt.Power(t.Lb_ISOLATION_PY,2)+mt.Power(t.Lb_ISOLATION_PZ,2)+mt.Power(m1old,2))
                    E2old = mt.Sqrt(mt.Power(t.Lb_ISOLATION_PX2,2)+mt.Power(t.Lb_ISOLATION_PY2,2)+mt.Power(t.Lb_ISOLATION_PZ2,2)+mt.Power(m2old,2))

                    Emu = mt.Sqrt(mt.Power(t.mu_PX,2)+mt.Power(t.mu_PY,2)+mt.Power(t.mu_PZ,2)+m_mu*m_mu)
                    
                    ETOT = E1 + E2 + Emu + ELc
                    ETOTold = E1old + E2old + Emu + ELc
                    pTOT_x = pLc12_x+t.mu_PX
                    pTOT_y = pLc12_y+t.mu_PY
                    pTOT_z = pLc12_z+t.mu_PZ
                    mTOT = mt.Sqrt(mt.Power(ETOT,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)
                    mTOTold = mt.Sqrt(mt.Power(ETOTold,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)



    
                #I save only those events with mLc12>2700 and with mTOT<5620
                    if mLc12>2700 and mTOT<5620 and (m1==m_K or m2==m_K):
                        isKenriched_new=1

            
                    h_E1pi.Fill(E1pi)
                    h_E2pi.Fill(E2pi)
                    h_ELc.Fill(ELc)
                    h_mLc12.Fill(mLc12)
                    h_mTOT_old.Fill(mTOTold)
                    h_mTOT_new.Fill(mTOT)
                    h_m1.Fill(m1)
                    h_m2.Fill(m2)
            if isKenriched_new==1:
                h_new.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,t.Event_PIDCalibEffWeight)

    print('Number of events selected with standard Kenriched selection: ',h_std.Integral())
    print('Number of events selected with new Kenriched selection: ',h_new.Integral())
    print('Number of entries with new Kenriched selection: ',h_new.GetEntries())

    h_std = ScaleHisto(h_std,1)
    if h_new.Integral()!=0:
        h_new = ScaleHisto(h_new,1)

        h_std.SetLineColor(r.kBlack)
        h_new.SetLineColor(r.kRed)

        #print('% of true kaons with PIDK>4.: ',h_PIDK_1.Integral()/h_PIDK.Integral())
        #print('% of true kaons with PIDK2>4.: ',h_PIDK2_1.Integral()/h_PIDK2.Integral())

        c = r.TCanvas('c','c',1500,500)
        c.Divide(3,1)
        c.cd(1)
        ymax = max([h_std.ProjectionX().GetBinContent(h_std.ProjectionX().GetMaximumBin()),h_new.ProjectionX().GetBinContent(h_new.ProjectionX().GetMaximumBin())])
        h_std.ProjectionX().GetYaxis().SetRangeUser(0,ymax+0.1*ymax)
        h_std.ProjectionX().Draw('hist')
        h_new.ProjectionX().Draw('hist same')
        c.cd(2)
        ymax = max([h_std.ProjectionY().GetBinContent(h_std.ProjectionY().GetMaximumBin()),h_new.ProjectionY().GetBinContent(h_new.ProjectionY().GetMaximumBin())])
        h_std.ProjectionY().GetYaxis().SetRangeUser(0,ymax+0.1*ymax)
        h_std.ProjectionY().Draw('hist')
        h_new.ProjectionY().Draw('hist same')
        c.cd(3)
        ymax = max([h_std.ProjectionZ().GetBinContent(h_std.ProjectionZ().GetMaximumBin()),h_new.ProjectionZ().GetBinContent(h_new.ProjectionZ().GetMaximumBin())])
        h_std.ProjectionZ().GetYaxis().SetRangeUser(0,ymax+0.1*ymax)
        h_std.ProjectionZ().Draw('hist')
        h_new.ProjectionZ().Draw('hist same')
        l = r.TLegend(0.1,0.7,0.4,0.9)
        l.AddEntry(h_std,'Standard','l')
        l.AddEntry(h_new,'New','l')
        l.Draw()
        c.SaveAs('plots/ComparisonKenrichedSel_'+sample+'.pdf')

        c1 = r.TCanvas('c1','c1',1000,500)
        c1.Divide(2,1)
        c1.cd(1)
        h_TrueID.Draw('hist')
        c1.cd(2)
        h_TrueID_1.Draw('hist')
        c1.SaveAs('plots/TrueIDparticle1_'+sample+'.pdf')

        c2 = r.TCanvas('c2','c2',1000,500)
        c2.Divide(2,1)
        c2.cd(1)
        h_TrueID2.Draw('hist')
        c2.cd(2)
        h_TrueID2_1.Draw('hist')
        c2.SaveAs('plots/TrueIDparticle2_'+sample+'.pdf')

        c3 = r.TCanvas('c3','c3',1000,500)
        c3.Divide(2,1)
        c3.cd(1)
        h_PIDK.Draw('hist')
        h_PIDK_1.SetFillColor(r.kRed)
        h_PIDK_1.SetFillStyle(3004)
        h_PIDK_1.Draw('hist same')
        c3.cd(2)
        h_PIDK2.Draw('hist')
        h_PIDK2_1.SetFillColor(r.kRed)
        h_PIDK2_1.SetFillStyle(3004)
        h_PIDK2_1.Draw('hist same')
        c3.SaveAs('plots/PIDK_'+sample+'.pdf')

        c9 = r.TCanvas('c9','c9',1500,500)
        c9.Divide(3,1)
        c9.cd(1)
        h_pK.SetLineColor(r.kRed)
        h_pK_1.SetLineColor(r.kBlack)
        h_pK.ProjectionX().Draw('hist')
        h_pK_1.ProjectionX().Draw('hist sames')
        c9.cd(2)
        h_pK.ProjectionY().Draw('hist')
        h_pK_1.ProjectionY().Draw('hist sames')
        c9.cd(3)
        h_pK.ProjectionZ().Draw('hist')
        h_pK_1.ProjectionZ().Draw('hist sames')
        l1 = r.TLegend(0.5,0.7,0.8,0.9)
        l1.AddEntry(h_pK,'True K ID','l')
        l1.AddEntry(h_pK_1,'PIDK>4','l')
        l1.Draw()
        c9.SaveAs('plots/MomentumK1_'+sample+'.pdf')


        c10 = r.TCanvas('c10','c10',1500,500)
        c10.Divide(3,1)
        c10.cd(1)
        h_pK2.SetLineColor(r.kRed)
        h_pK2_1.SetLineColor(r.kBlack)
        h_pK2.ProjectionX().Draw('hist')
        h_pK2_1.ProjectionX().Draw('hist sames')
        c10.cd(2)
        h_pK2.ProjectionY().Draw('hist')
        h_pK2_1.ProjectionY().Draw('hist sames')
        c10.cd(3)
        h_pK2.ProjectionZ().Draw('hist')
        h_pK2_1.ProjectionZ().Draw('hist sames')
        l2 = r.TLegend(0.5,0.7,0.8,0.9)
        l2.AddEntry(h_pK2,'True K ID (2)','l')
        l2.AddEntry(h_pK2_1,'PIDK2>4','l')
        l2.Draw()
        c10.SaveAs('plots/MomentumK2_'+sample+'.pdf')

        c11= r.TCanvas('c11','c11',1000,500)
        c11.Divide(2,1)
        c11.cd(1)
        h_PIDKp.Draw('colz')
        c11.cd(2)
        h_PIDK2p2.Draw('colz')
        c11.SaveAs('plots/PIDvsMom_'+sample+'.pdf')
        
        c8 = r.TCanvas('c8','c8',1000,500)
        c8.Divide(2,1)
        c8.cd(1)
        h_PIDdiff.Draw()
        h_PIDdiff_1.SetFillColor(r.kRed)
        h_PIDdiff_1.SetFillStyle(3004)
        h_PIDdiff_1.Draw('sames')
        c8.cd(2)
        h_PIDdiff2.Draw()
        h_PIDdiff2_1.SetFillColor(r.kRed)
        h_PIDdiff2_1.SetFillStyle(3004)
        h_PIDdiff2_1.Draw('sames')
        c8.SaveAs('plots/PIDp-PIDK_'+sample+'.pdf')

        c13 = r.TCanvas('c13','c13',500,500)
        h_mTOT_new.SetLineColor(r.kRed)
        h_mTOT_new.Draw()
        h_mTOT_old.SetLineColor(r.kBlack)
        h_mTOT_old.Draw('sames')
        l = r.TLegend(0.1,0.7,0.4,0.9)
        l.AddEntry(h_mTOT_old,'Standard','l')
        l.AddEntry(h_mTOT_new,'New','l')
        l.Draw()
        c13.SaveAs('plots/RecoLbmTOT_'+sample+'.pdf')


    '''
    c4 = r.TCanvas('c4','c4',500,500)
    h_mLc12.Draw()

    c5 = r.TCanvas('c5','c5',1500,500)
    c5.Divide(3,1)
    c5.cd(1)
    h_E1pi.Draw('hist')
    c5.cd(2)
    h_E2pi.Draw('hist')
    c5.cd(3)
    h_ELc.Draw('hist')
    '''

    c6 = r.TCanvas('c6','c6',1000,500)
    c6.Divide(2,1)
    c6.cd(1)
    #h_m1.Draw('hist')
    t.Draw('Lb_ISOLATION_PIDK')
    c6.cd(2)
    t.Draw('Lb_ISOLATION_PIDK2')
    #h_m2.Draw('hist')

'''
    c7 = r.TCanvas('c7','c7',500,500)
    h_muCh.Draw()

    '''
