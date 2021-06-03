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
picodes = [111,211]
pcodes = [2212]

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h


#Upload files to read
datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
#samples = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
#samples = ['Lctaunu','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
#samples = ['Lcmunu','LcDs']
samples = ['LcDs']

blist = ['FinalSel','isKenriched','Lb_ISOLATION_TruePID','Lb_ISOLATION_TrueMotherPID','Lb_ISOLATION_TruePID2','Lb_ISOLATION_TrueMotherPID2','Lb_ISOLATION_BDT','Lb_ISOLATION_BDT2','Lb_ISOLATION_BDT3','Lb_ISOLATION_PIDK','Lb_ISOLATION_PIDp','Lb_ISOLATION_PIDK2','Lb_ISOLATION_PIDp2','Lb_ISOLATION_PIDK3','Lb_ISOLATION_PIDp3','Lb_ISOLATION_Type', 'Lb_ISOLATION_Type2','Lb_ISOLATION_Type3','Lb_ISOLATION_NNghost', 'Lb_ISOLATION_NNghost2','Lb_ISOLATION_NNghost3', 'Lb_ISOLATION_CHARGE', 'Lb_ISOLATION_CHARGE2','Lb_ISOLATION_CHARGE3','Event_PIDCalibEffWeight','mu_ID','Lb_ISOLATION_PX','Lb_ISOLATION_PX2','Lb_ISOLATION_PY','Lb_ISOLATION_PY2','Lb_ISOLATION_PZ','Lb_ISOLATION_PZ2','Lc_PX','Lc_PY','Lc_PZ','mu_PX','mu_PY','mu_PZ','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc']

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
    #and turn on only those needed
    for b in blist:
        t.SetBranchStatus(b,1)

    #counting variables
    nsamech, nsamech2=0,0 #number of anti-iso particles with charge = to that of the muon
    noppch, noppch2=0,0 #number of anti-iso particles with charge opposite to that of the muon
    nneutr,nneutr2=0,0 #number of anti-iso particles with 0 charge
    nPIDKg4, nPIDK2g4=0,0 #number of anti-iso particles with PIDK>4
    nssPIDKg4, nssPIDK2g4=0,0 #number of anti-iso particles with PIDK>4 (same ch)
    nosPIDKg4, nosPIDK2g4=0,0 #number of anti-iso particles with PIDK>4 (opp ch)
    nssPIDKg4PIDp, nssPIDK2g4PIDp=0,0 #number of anti-iso particles with PIDK>4 (same ch) and PIDp-PIDK<=0
    nClassK, nClassK2 = 0,0 #number of anti-iso particles classified as Kaons
    nClassp, nClassp2 = 0,0 #number of anti-iso particles classified as protons
    nClasspi, nClasspi2 = 0,0 #number of anti-iso particles classified as pions
    
    nKPIDKg4, nKPIDK2g4=0,0 #number of anti-iso true K with PIDK>4
    nTrueK, nTrueK2=0,0 #number of True K
    nTruep, nTruep2=0,0 #number of True protons
    nTruepi, nTruepi2=0,0 #number of True pions
    nOther, nOther2 =0,0 #number of anti-iso particles which are nor pi, p or K
    nKss, nK2ss, nKos, nK2os = 0,0,0,0 #number of true K with same (opposite) sign wrt to mu
    nKssPIDKg4, nK2ssPIDK2g4=0,0 #number of anti-iso true K (same ch) with PIDK>4
    nKssPIDKg4PIDp, nK2ssPIDK2g4PIDp=0,0 #number of anti-iso true K (same ch) with PIDK>4 and PIDp-PIDk<=0
    nKosPIDKg4, nK2osPIDK2g4=0,0 #number of anti-iso true K (opposite ch) with PIDK>4

    h_KMID = r.TH1F('h_KMID','Mother ID True Kaons;TrueID;',10000,0,5000)
    h_K2MID = r.TH1F('h_K2MID','Mother ID True Kaons (2nd particle);TrueID;',10000,0,5000)

    #Loop over entries
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if t.isIsolated==1 or t.FinalSel!=1:
            continue
        muCh = -t.mu_ID/13
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

        if BDT>ISOBDTcut and BDT2>ISOBDT2cut:
            if Type==3 and NNghost<0.2 and Type2==3 and NNghost2<0.2:
                if Ch==muCh:
                    nsamech+=1
                if Ch==-muCh:
                    noppch+=1
                if Ch2==muCh:
                    nsamech2+=1
                if Ch2==-muCh:
                    noppch2+=1
                if PIDK>4:
                    nPIDKg4+=1
                    if Ch==muCh:
                        nssPIDKg4+=1
                        if PIDp-PIDK<=0:
                            nssPIDKg4PIDp+=1
                        else:
                            nClassp+=1
                    if Ch==-muCh:
                        nosPIDKg4+=1
                else:
                    nClasspi+=1

                if PIDK2>4:
                    nPIDK2g4+=1
                    if Ch2==muCh:
                        nssPIDK2g4+=1
                        if PIDp2-PIDK2<=0:
                            nssPIDK2g4PIDp+=1
                        else:
                            nClassp2+=1
                    if Ch2==-muCh:
                        nosPIDK2g4+=1
                else:
                    nClasspi2+=1
                
                if r.TMath.Abs(TrueID) in Kcodes:
                    nTrueK+=1
                    h_KMID.Fill(r.TMath.Abs(TrueMID))
                    if Ch==muCh:
                        nKss+=1
                        if PIDK>4:
                            nKssPIDKg4+=1
                            if PIDp-PIDK<=0:
                                nKssPIDKg4PIDp+=1
                    if Ch==-muCh:
                        nKos+=1
                        if PIDK>4:
                            nKosPIDKg4+=1

                if r.TMath.Abs(TrueID2) in Kcodes:
                    h_K2MID.Fill(r.TMath.Abs(TrueMID2))
                    nTrueK2+=1
                    if Ch2==muCh:
                        nK2ss+=1
                        if PIDK2>4:
                            nK2ssPIDK2g4+=1
                            if PIDp2-PIDK2<=0:
                                nK2ssPIDK2g4PIDp+=1
                    if Ch2==-muCh:
                        nK2os+=1
                        if PIDK2>4:
                            nK2osPIDK2g4+=1
                
                if r.TMath.Abs(TrueID) in pcodes:
                    nTruep+=1
                
                if r.TMath.Abs(TrueID) in picodes:
                    nTruepi+=1
                
                if r.TMath.Abs(TrueID2) in pcodes:
                    nTruep2+=1
                
                if r.TMath.Abs(TrueID2) in picodes:
                    nTruepi2+=1


    nClassK = nssPIDKg4PIDp+ nosPIDKg4
    nClassK2 = nssPIDK2g4PIDp+ nosPIDK2g4
    nKPIDKg4 = nKosPIDKg4+nKssPIDKg4 
    nKPIDK2g4 = nK2osPIDK2g4+nK2ssPIDK2g4 

    print('Number of anti-isolated particles with same charge as the muon (1 and 2): ',nsamech, nsamech2)
    print('Number of anti-isolated particles with opposite charge to the muon (1 and 2): ',noppch, noppch2)
    print('Number of neutral anti-isolated particles (1 and 2): ',nneutr, nneutr2)
    print('-------------------------------------------------------------')
    print('              Values estimated using true info')
    print('-------------------------------------------------------------')
    print('Number of anti-isolated true K: ', nTrueK, nTrueK2)
    print('Number of anti-isolated true pi: ', nTruepi, nTruepi2)
    print('Number of anti-isolated true protons: ', nTruep, nTruep2)
    print('Number of anti-isolated true K with same charge of the muon: ', nKss, nK2ss)
    print('Number of anti-isolated true K with opposite charge of the muon: ', nKos, nK2os)
    print('Number of anti-isolated true K with PIDK>4: ', nKPIDKg4, nKPIDK2g4)
    print('Number of anti-isolated true K with opposite ch and PIDK>4: ', nKosPIDKg4, nK2osPIDK2g4)
    print('Number of anti-isolated true K with same ch and PIDK>4: ', nKssPIDKg4, nK2ssPIDK2g4)
    print('Number of anti-isolated true K with same ch and PIDK>4 and PIDp-PIDK<=0: ', nKssPIDKg4PIDp, nK2ssPIDK2g4PIDp)
    print('-------------------------------------------------------------')
    print('              Values estimated using PID only')
    print('-------------------------------------------------------------')
    print('Number of anti-isolated particles with PIDK>4: ', nPIDKg4, nPIDK2g4)
    print('Number of anti-isolated particles with PIDK>4 (opp ch): ', nosPIDKg4, nosPIDK2g4)
    print('Number of anti-isolated particles with PIDK>4 (same ch): ', nssPIDKg4, nssPIDK2g4)
    print('Number of anti-isolated particles with PIDK>4 and PIDp-PIDK<=0 (same ch): ', nssPIDKg4PIDp, nssPIDK2g4PIDp)
    print('Number of anti-isolated particles classified as K: ', nClassK, nClassK2)
    print('Number of anti-isolated particles classified as pions: ', nClasspi, nClasspi2)
    print('Number of anti-isolated particles classified as protons: ', nClassp, nClassp2)

    c = r.TCanvas()
    c.Divide(2,1)
    c.cd(1)
    h_KMID.Draw()
    c.cd(2)
    h_K2MID.Draw()


