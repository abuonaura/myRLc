import ROOT as r

files = {'MagUp':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root','MagDown':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root'}
csamples = {'MagUp':'$FILEDIR/ControlSamples/MC_Kenr/Lb_LcDs_MagUp.root','MagDown':'$FILEDIR/ControlSamples/MC_Kenr/Lb_LcDs_MagDown.root'}
polarities=['MagUp','MagDown']

def GetTotalNumberMCevts():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.GetEntries()
    print('Total number of MC events before applying ANY selection: ', ntotMC)
    return ntotMC

def GetTotalNumber2body():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lb_TrueHadron_D2_ID==0','goff')
        '''
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  t.Lb_TrueHadron_D2_ID==0:
                ntot_2body+=1
        '''
    print('Total number of 2body MC events before applying ANY selection: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lb_TrueHadron_D2_ID!=0','goff')
    print('Total number of mbody MC events before applying ANY selection: ', ntot_mbody)
    return ntot_mbody

def GetTotalNumberMCevts_TMatch():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50','goff')
    print('Total number of MC events after applying Truth Matching: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID==0','goff')
    print('Total number of 2body MC events after applying Truth Matching: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID!=0','goff')
    print('Total number of mbody MC events after applying Truth Matching: ', ntot_mbody)
    return ntot_mbody

def GetTotalNumberMCevts_TMatch_Iso():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of MC events after applying Truth Matching + Isolation: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch_Iso():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID==0&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of 2body MC events after applying Truth Matching + Isolation: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch_Iso():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID!=0&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of mbody MC events after applying Truth Matching + Isolation: ', ntot_mbody)
    return ntot_mbody

def GetTotalNumberMCevts_TMatch_Kenriched():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50','goff')
    print('Total number of MC events after applying Truth Matching +  Kenriched selection: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch_Kenriched():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID==0','goff')
    print('Total number of 2body MC events after applying Truth Matching +  Kenriched selection: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch_Kenriched():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID!=0','goff')
    print('Total number of mbody MC events after applying Truth Matching + Kenriched selection : ', ntot_mbody)
    return ntot_mbody

def FillHistoD1(h,pdg):
    if r.TMath.Abs(pdg)==411:
        h.Fill(1)
    if r.TMath.Abs(pdg)==421:
        h.Fill(2)
    if r.TMath.Abs(pdg)==413:
        h.Fill(3)
    if r.TMath.Abs(pdg)==423:
        h.Fill(4)
    if r.TMath.Abs(pdg)==415:
        h.Fill(5)
    if r.TMath.Abs(pdg)==425:
        h.Fill(6)
    if r.TMath.Abs(pdg)==431:
        h.Fill(7)
    if r.TMath.Abs(pdg)==10431:
        h.Fill(8)
    if r.TMath.Abs(pdg)==433:
        h.Fill(9)
    if r.TMath.Abs(pdg)==435:
        h.Fill(10)
    h.SetDirectory(0)
    return h

def FillHistoD2(h,pdg):
    if r.TMath.Abs(pdg)==311:
        h.Fill(1)
    if r.TMath.Abs(pdg)==321:
        h.Fill(2)
    if r.TMath.Abs(pdg)==10311:
        h.Fill(3)
    h.SetDirectory(0)
    return h



def Check_mbody():
    h_BDT = r.TH1F('h_BDT','Lb_ISOLATION_BDT',50,-2,1)
    h_BDT_sel = r.TH1F('h_BDT_sel','Lb_ISOLATION_BDT',50,-2,1)
    h_BDT2 = r.TH1F('h_BDT2','Lb_ISOLATION_BDT2',50,-2,1)
    h_BDT2_sel = r.TH1F('h_BDT2_sel','Lb_ISOLATION_BDT2',50,-2,1)
    h_TrueD1 = r.TH1F('h_TrueD1','Lb_TrueHadron_D1_ID',10,1,11)
    h_TrueD1_1= r.TH1F('h_TrueD1_1','Lb_TrueHadron_D1_ID',10,1,11) #PASSING ISO_BDT>0.35 AND ISO_BDT2>0.2
    D1particles = ['D^{+}','D^{0}','D^{*}(2010)^{+}','D^{*}(2007)^{0}','D^{*}_{2}(2460)^{+}','D^{*}_{2}(2460)^{0}','D^{+}_{s}','D^{*}_{s0}(2317)^{+}','D^{*+}_{s}','D^{*}_{s}(2573)^{+}']
    for n, d1p in enumerate(D1particles):
        print(n, d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1_1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1_1.GetXaxis().LabelsOption('v')
        h_TrueD1_1.GetXaxis().LabelsOption('v')
    h_TrueD2 = r.TH1F('h_TrueD2','Lb_TrueHadron_D2_ID',3,1,4)
    h_TrueD2_1 = r.TH1F('h_TrueD2_1','Lb_TrueHadron_D2_ID',3,1,4)  #PASSING ISO_BDT>0.35 AND ISO_BDT2>0.2
    D2particles = ['K^{0}','K','K^{0*}']
    for n, d2p in enumerate(D2particles):
        print(n, d2p)
        h_TrueD2.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2.GetXaxis().LabelsOption('v')
        h_TrueD2_1.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2_1.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2_1.GetXaxis().LabelsOption('v')
    h_PIDK = r.TH1F('h_PIDK','Lb_ISOLATION_PIDK',50,-200,200)
    h_PIDK_1 = r.TH1F('h_PIDK_1','Lb_ISOLATION_PIDK',50,-200,200) #PASSING BOTH ISO BDT CUTS 
    h_PIDK_2 = r.TH1F('h_PIDK_2','Lb_ISOLATION_PIDK',50,-200,200) #PASSING BOTH ISO BDT CUTS + PIDK>4 (= K)
    h_PIDK2 = r.TH1F('h_PIDK2','Lb_ISOLATION_PIDK2',50,-200,200)
    h_PIDK2_1 = r.TH1F('h_PIDK2_1','Lb_ISOLATION_PIDK2',50,-200,200) #PASSING BOTH ISO BDT CUTS 
    h_PIDK2_2 = r.TH1F('h_PIDK2_2','Lb_ISOLATION_PIDK2',50,-200,200) #PASSING BOTH ISO BDT CUTS + PIDK>4 (= K)

    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  t.Lb_TrueHadron_D2_ID!=0 and t.Lc_BKGCAT<30 and t.Lb_BKGCAT<50:
                h_BDT.Fill(t.Lb_ISOLATION_BDT)
                h_BDT2.Fill(t.Lb_ISOLATION_BDT2)
                h_TrueD1=FillHistoD1(h_TrueD1,t.Lb_TrueHadron_D1_ID)
                h_TrueD2=FillHistoD2(h_TrueD2,t.Lb_TrueHadron_D2_ID)
                h_PIDK.Fill(t.Lb_ISOLATION_PIDK)
                h_PIDK2.Fill(t.Lb_ISOLATION_PIDK2)
                if t.Lb_ISOLATION_BDT>0.35:
                    h_BDT_sel.Fill(t.Lb_ISOLATION_BDT)
                if t.Lb_ISOLATION_BDT2>0.2:
                    h_BDT2_sel.Fill(t.Lb_ISOLATION_BDT2)
                if t.Lb_ISOLATION_BDT>0.35 and t.Lb_ISOLATION_BDT2>0.2:
                    h_TrueD1_1=FillHistoD1(h_TrueD1_1,t.Lb_TrueHadron_D1_ID)
                    h_TrueD2_1=FillHistoD2(h_TrueD2_1,t.Lb_TrueHadron_D2_ID)
                    h_PIDK_1.Fill(t.Lb_ISOLATION_PIDK)
                    h_PIDK2_1.Fill(t.Lb_ISOLATION_PIDK2)
                    if t.Lb_ISOLATION_PIDK>4:
                        h_PIDK_2.Fill(t.Lb_ISOLATION_PIDK)
                    if t.Lb_ISOLATION_PIDK2>4:
                        h_PIDK2_2.Fill(t.Lb_ISOLATION_PIDK2)


    c = r.TCanvas('c','c',1500,1500)
    c.Divide(4,2)
    c.cd(1)
    h_BDT.Draw()
    h_BDT_sel.SetFillStyle(3004)
    h_BDT_sel.SetFillColor(r.kRed)
    h_BDT_sel.Draw('same')
    print('Fraction of events with ISO_BDT > 0.35: ', h_BDT_sel.GetEntries()*1./h_BDT.GetEntries())
    c.cd(2)
    h_BDT2.Draw()
    h_BDT2_sel.SetFillStyle(3004)
    h_BDT2_sel.SetFillColor(r.kRed)
    h_BDT2_sel.Draw('same')
    print('Fraction of events with ISO_BDT2 > 0.2: ', h_BDT2_sel.GetEntries()*1./h_BDT2.GetEntries())
    c.cd(3)
    h_TrueD1.Draw()
    h_TrueD1_1.SetFillStyle(3004)
    h_TrueD1_1.SetFillColor(r.kRed)
    h_TrueD1_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2: ',  h_TrueD1_1.GetEntries()*1/h_TrueD1.GetEntries())
    c.cd(4)
    h_TrueD2.Draw()
    h_TrueD2_1.SetFillStyle(3004)
    h_TrueD2_1.SetFillColor(r.kRed)
    h_TrueD2_1.Draw('same')
    c.cd(5)
    h_PIDK.Draw()
    h_PIDK_1.SetFillStyle(3004)
    h_PIDK_1.SetFillColor(r.kRed)
    h_PIDK_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2 and PIDK>4: ', h_PIDK_2.GetEntries()*1./ h_PIDK_1.GetEntries())
    c.cd(6)
    h_PIDK2.Draw()
    h_PIDK2_1.SetFillStyle(3004)
    h_PIDK2_1.SetFillColor(r.kRed)
    h_PIDK2_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2 and PIDK2>4: ', h_PIDK2_2.GetEntries()*1./ h_PIDK2_1.GetEntries())
    c.Update()
    c.Draw()

    return 


GetTotalNumberMCevts()
GetTotalNumber2body()
GetTotalNumberMbody()
print()
GetTotalNumberMCevts_TMatch()
GetTotalNumber2body_TMatch()
GetTotalNumberMbody_TMatch()
print()
GetTotalNumberMCevts_TMatch_Iso()
GetTotalNumber2body_TMatch_Iso()
GetTotalNumberMbody_TMatch_Iso()
print()
GetTotalNumberMCevts_TMatch_Kenriched()
GetTotalNumber2body_TMatch_Kenriched()
GetTotalNumberMbody_TMatch_Kenriched()
print()
Check_mbody()
