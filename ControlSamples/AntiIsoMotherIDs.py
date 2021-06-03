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

    #Loop over entries
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if t.FinalSel!=1:
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
             print(TrueMID, TrueMID2)
