import ROOT as r

datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'

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


