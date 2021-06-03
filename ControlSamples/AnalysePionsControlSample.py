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
PiMotherIDs = {'Lcmunu':4122,'Lctaunu':4122,'LcDs':4122,'Lc2625munu':104124,'Lc2625taunu':104124,'Lc2625Ds':104124,'Lc2593munu':14122,'Lc2593taunu':14122,'Lc2593Ds':14122}


blist = ['FinalSel','isKenriched','Lb_ISOLATION_TruePID','Lb_ISOLATION_TrueMotherPID','Lb_ISOLATION_TruePID2','Lb_ISOLATION_TrueMotherPID2','Lb_ISOLATION_BDT','Lb_ISOLATION_BDT2','Lb_ISOLATION_BDT3','Lb_ISOLATION_PIDK','Lb_ISOLATION_PIDp','Lb_ISOLATION_PIDK2','Lb_ISOLATION_PIDp2','Lb_ISOLATION_PIDK3','Lb_ISOLATION_PIDp3','Lb_ISOLATION_Type', 'Lb_ISOLATION_Type2','Lb_ISOLATION_Type3','Lb_ISOLATION_NNghost', 'Lb_ISOLATION_NNghost2','Lb_ISOLATION_NNghost3', 'Lb_ISOLATION_CHARGE', 'Lb_ISOLATION_CHARGE2','Lb_ISOLATION_CHARGE3','Event_PIDCalibEffWeight','mu_ID','Lb_ISOLATION_PX','Lb_ISOLATION_PX2','Lb_ISOLATION_PY','Lb_ISOLATION_PY2','Lb_ISOLATION_PZ','Lb_ISOLATION_PZ2','Lc_PX','Lc_PY','Lc_PZ','mu_PX','mu_PY','mu_PZ','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc']


def FillMotherHisto(h,pdg,w):
    if r.TMath.Abs(pdg)==113: #rho0
        h.Fill(1,w)
    elif r.TMath.Abs(pdg)==213: #rho+
        h.Fill(2,w)
    elif r.TMath.Abs(pdg)==221: #eta
        h.Fill(3,w)
    elif r.TMath.Abs(pdg)==331: #eta'
        h.Fill(4,w)
    elif r.TMath.Abs(pdg)==223: #omega(782)
        h.Fill(5,w)
    elif r.TMath.Abs(pdg)==333: #phi
        h.Fill(6,w)
    elif r.TMath.Abs(pdg)==310: #K0S
        h.Fill(7,w)
    elif r.TMath.Abs(pdg)==411: #D+
        h.Fill(8,w)
    elif r.TMath.Abs(pdg)==421: #D0
        h.Fill(9,w)
    elif r.TMath.Abs(pdg)==413: #D*+
        h.Fill(10,w)
    elif r.TMath.Abs(pdg)==423: #D*0
        h.Fill(11,w)
    elif r.TMath.Abs(pdg)==431: #Ds
        h.Fill(12,w)
    elif r.TMath.Abs(pdg)==433: #D*s
        h.Fill(13,w)
    elif r.TMath.Abs(pdg)==4122:#Lc
        h.Fill(14,w)
    elif r.TMath.Abs(pdg)==104124:#Lc(2625)
        h.Fill(15,w)
    elif r.TMath.Abs(pdg)==14122: #Lc(2593)
        h.Fill(16,w)
    else:
        h.Fill(17,w)
    h.SetDirectory(0)
    return h

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

    hBDT = r.TH1F('hBDT','ISOBDT',25,-1,1) #Iso BDT pi from LcStar
    hBDTo = r.TH1F('hBDTo','ISOBDT',25,-1,1)  #Iso BDT pi NOT from LcStar
    hBDTall = r.TH1F('hBDTall','ISOBDT',25,-1,1) #ISOBDT all pions
    hBDT2 = r.TH1F('hBDT2','ISOBDT2',25,-1,1)
    hBDT2o = r.TH1F('hBDT2o','ISOBDT2',25,-1,1)
    hBDT2all = r.TH1F('hBDT2all','ISOBDT2',25,-1,1)

    nTruePiPi = 0
    nPiPiall=0

    hMotherID = r.TH1F('hMotherID','MotherID True Pi (AntiIso 1)',17,1,18)
    hMotherID2 = r.TH1F('hMotherID2','MotherID True Pi (AntiIso 2)',17,1,18)
    MotherParticles = ['#rho^{0}','#rho^+','#eta','#eta^{1}','#omega(782)','#phi','K^{0}_{S}','D^{+}','D^{0}','D^{*+}','D^{*0}','D_{s}','D^{*}_{s}','#Lambda_{c}','#Lambda_{c}(2625)','#Lambda_{c}(2593)','others']
    for n,mp in enumerate(MotherParticles):
        hMotherID.GetXaxis().SetBinLabel(n+1,mp)
        hMotherID.GetXaxis().LabelsOption('v')
        hMotherID2.GetXaxis().SetBinLabel(n+1,mp)
        hMotherID2.GetXaxis().LabelsOption('v')




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

        if Type==3 and NNghost<0.2 and Type2==3 and NNghost2<0.2:
            if abs(TrueID)==211 and abs(TrueID2)==211:
                nPiPiall +=t.Event_PIDCalibEffWeight
                hMotherID = FillMotherHisto(hMotherID,TrueMID,t.Event_PIDCalibEffWeight)
                hMotherID2 = FillMotherHisto(hMotherID2,TrueMID2,t.Event_PIDCalibEffWeight)
                hBDTall.Fill(BDT,t.Event_PIDCalibEffWeight)
                hBDT2all.Fill(BDT2,t.Event_PIDCalibEffWeight)
                if abs(TrueMID)==PiMotherIDs[sample] and abs(TrueMID2)==PiMotherIDs[sample]:
                    nTruePiPi+=t.Event_PIDCalibEffWeight
                    hBDT.Fill(BDT,t.Event_PIDCalibEffWeight)
                    hBDT2.Fill(BDT2,t.Event_PIDCalibEffWeight)
                else:
                    hBDTo.Fill(BDT,t.Event_PIDCalibEffWeight)
                    hBDT2o.Fill(BDT2,t.Event_PIDCalibEffWeight)

    c = r.TCanvas('c','',1000,500)
    c.Divide(2,1)
    c.cd(1)
    hBDTall.SetMarkerColor(r.kBlack)
    hBDTall.SetMarkerStyle(20)
    hBDTall.Draw('E0')
    hBDTo.SetLineColor(r.kRed)
    hBDTo.Draw('same')
    hBDT.Draw('same')
    c.cd(2)
    hBDT2all.SetMarkerColor(r.kBlack)
    hBDT2all.SetMarkerStyle(20)
    hBDT2all.Draw('E0')
    hBDT2o.SetLineColor(r.kRed)
    hBDT2o.Draw('same')
    hBDT2.Draw('same')

    c1 = r.TCanvas('c1','MotherIDs',1000,500)
    c1.Divide(2,1)
    c1.cd(1)
    hMotherID.Draw()
    c1.cd(2)
    hMotherID2.Draw()
