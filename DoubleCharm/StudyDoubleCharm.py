import ROOT as r

datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'

fname = datadir+'Lb_LcDs_MagUp_full.root'
preselfname = datadir+'Lb_LcDs_MagUp_full_preselectionVars.root'
pidcalibFile = fname[0:-5]+'_PIDCalib.root'

f = r.TFile(fname,'READ')
t = f.Get('tupleout/DecayTree')
fpresel = r.TFile(preselfname,'READ')
tpresel = fpresel.Get('DecayTree')
t.AddFriend(tpresel)
fpidcalib = r.TFile(pidcalibFile)
tcalib = fpidcalib.Get('tupleout/DecayTree')
t.AddFriend(tcalib)

h_ISOBDT = r.TH1F('h_ISOBDT','',50,-1,1)
h_ISOBDT2 = r.TH1F('h_ISOBDT2','',50,-1,1)
h_ISOBDT_2 = r.TH1F('h_ISOBDT_2','',50,-1,1)
h_ISOBDT2_2 = r.TH1F('h_ISOBDT2_2','',50,-1,1)

m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_e = 0.548 #MeV

def PlotIsoBDTForDifferentDecays():
    for i in range(500000):
        t.GetEntry(i)
        if t.FinalSel==True:
            if i%100000==0:
                print(i)
            weight = t.Event_PIDCalibEffWeight
            LbD1 = t.Lb_TrueHadron_D1_ID #2nd daughter Lb (1st is Lc)
            LbD2 = t.Lb_TrueHadron_D2_ID
            LbD1GD0 = t.Lb_TrueHadron_D1_GD0_ID
            LbD1GD1 = t.Lb_TrueHadron_D1_GD1_ID
            LbD1GD2 = t.Lb_TrueHadron_D1_GD2_ID
            LbD2GD0 = t.Lb_TrueHadron_D2_GD0_ID
            LbD2GD1 = t.Lb_TrueHadron_D2_GD1_ID
            LbD2GD2 = t.Lb_TrueHadron_D2_GD2_ID
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
            #I want to investigate decay in LcD0(->K/K*/pi/mu)K+
            if abs(LbD1)==421 and abs(LbD2)==321:
                if (abs(LbD1GD0)==321 or abs(LbD1GD0)==323 or abs(LbD1GD0)==211) and abs(LbD1GD1)==13:
                    h_ISOBDT.Fill(BDT,weight)
                    h_ISOBDT2.Fill(BDT2,weight)
            #I want to investigate decay in LcDs(->eta(->pipipi0)mu) or phi
            if abs(LbD1)==431:
                if (abs(LbD1GD0)==221 or abs(LbD1GD0)==333) and abs(LbD1GD1)==13:
                    h_ISOBDT_2.Fill(BDT,weight)
                    h_ISOBDT2_2.Fill(BDT2,weight)
                
    c = r.TCanvas('c','c',1000,500)
    c.Divide(2,1)
    c.cd(1)
    h_ISOBDT.Draw('hist')
    c.cd(2)
    h_ISOBDT2.Draw('hist')

    c1 = r.TCanvas('c1','c1',1000,500)
    c1.Divide(2,1)
    c1.cd(1)
    h_ISOBDT_2.Draw('hist')
    c1.cd(2)
    h_ISOBDT2_2.Draw('hist')

def AssignMass(PID):
    m=0
    if abs(PID)==211: m=m_pi
    if abs(PID)==11: m=m_e
    if abs(PID)==13: m=m_mu
    if abs(PID)==321: m = m_K
    if abs(PID)==2212: m= m_p
    return m

def CheckHowManyTimesCorrectMassAssigned():
    nevts=0
    ncorrect_1, ncorrect_2 = 0,0
    m1,m2=0,0
    nm10, nm20, nm1true0, nm2true0 = 0,0,0,0
    ne1, nmu1, ne2, nmu2,np1True, np2True, nK1True, nK2True,nPi1True, nPi2True=0,0,0,0,0,0,0,0,0,0
    np1, np2, nK1, nK2,nPi1, nPi2=0,0,0,0,0,0
    #for i in range(t.GetEntries()):
    for i in range(100000):
        t.GetEntry(i)
        if t.FinalSel==False:
            continue
        if i%100000==0:
            print(i)
        PIDtrue = t.Lb_ISOLATION_TruePID
        m1true=AssignMass(PIDtrue)
        PIDtrue2 = t.Lb_ISOLATION_TruePID2
        m2true=AssignMass(PIDtrue2)
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
        if t.Lb_ISOLATION_BDT>0.35 and t.Lb_ISOLATION_BDT2>0.2:
            if abs(PIDtrue)==2212: np1True+=1 
            if abs(PIDtrue2)==2212: np2True+=1
            if abs(PIDtrue)==321: nK1True+=1
            if abs(PIDtrue2)==321: nK2True+=1
            if abs(PIDtrue)==211: nPi1True+=1
            if abs(PIDtrue2)==211: nPi2True+=1
            if abs(PIDtrue)==11: ne1+=1
            if abs(PIDtrue2)==11: ne2+=1
            if abs(PIDtrue)==13: nmu1+=1
            if abs(PIDtrue2)==13: nmu2+=1

            if PIDK>4.:
                if (Ch==-muCh):
                    #I also assume this particle is a K
                    m1 = m_K
                    nK1+=1
                if Ch==muCh and (PIDp-PIDK)<0.:
                    #I also assume this particle is a K
                    m1 = m_K
                    nK1+=1
                if Ch==muCh and (PIDp-PIDK)>0.:
                    #In this case I assume it is a proton
                    m1 = m_p
                    np1+=1
            else: #if it does not satisfy the hypothesis for being a K, I assume it is a pi
                m1 = m_pi
                nPi1+=1
            #and I check which masses I can assume for the second
            if PIDK2>4.:
                if (Ch2==-muCh or (Ch2==muCh and (PIDp2 - PIDK2)<0.)):
                    #I also assume this particle is a K
                    m2 = m_K
                    nK2+=1
                if (Ch2==muCh and PIDp2-PIDK2>0.):
                    #In this case I assume it is a proton
                    m2 = m_p
                    np2+=1
            else: #if it does not satisfy the hypothesis for being a K, I assume it is a pi
                m2 = m_pi
                nPi2+=1
            nevts+=1
            if m1==m1true:
                ncorrect_1+=1
            if m2==m2true:
                ncorrect_2+=1

    if m1==0:
        nm10+=1
    if m2==0:
        nm20+=1
    if m1true==0:
        nm1true0+=1
        print(PIDtrue)
    if m2true==0:
        nm2true0+=1
    print(nevts, ncorrect_1, ncorrect_2)
    print('% of correct selection particle 1', ncorrect_1/nevts)
    print('% of correct selection particle 2', ncorrect_2/nevts)
    print('nevts with m1=0: ',nm10)
    print('nevts with m1true=0: ',nm1true0)
    print('nevts with m2=0: ',nm20)
    print('nevts with m2true=0: ',nm2true0)
    print(np1True, np2True, nK1True, nK2True,nPi1True, nPi2True, ne1, ne2, nmu1, nmu2)
    print(np1, np2, nK1, nK2, nPi1, nPi2)
