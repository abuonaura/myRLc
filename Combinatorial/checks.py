import ROOT as r

f_dataSS_magUp = r.TFile('/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/ControlSamples/Lb_DataSS_MagUp_reduced_preselected_Kenr_sw.root','READ')
t = f_dataSS_magUp.Get('DecayTree')
f_dataSS_magDown = r.TFile('/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/ControlSamples/Lb_DataSS_MagDown_reduced_preselected_Kenr_sw.root','READ')
t1 = f_dataSS_magDown.Get('DecayTree')

h_El = r.TH1F('h_El','',10,0,2600)
h_q2 = r.TH1F('h_q2','',4,-2E6,14E6)
h_Mmiss = r.TH1F('h_Mmiss2','',10,-2E6,14E6)
h_mLc = r.TH1F('h_mLc','',100,2230,2330)
h_dataSS = r.TH3F('h_dataSS','; q2^{2} (MeV^{2}); E_{l} (MeV^{2}); Mmiss^{2} (MeV^{2})',4,-2E6,14E6,10,0,2600,10,-2E6,14E6)
h1_El = r.TH1F('h1_El','',10,0,2600)
h1_q2 = r.TH1F('h1_q2','',4,-2E6,14E6)
h1_Mmiss = r.TH1F('h1_Mmiss2','',10,-2E6,14E6)
h1_mLc = r.TH1F('h1_mLc','',100,2230,2330)
h_Pi = r.TH3F('h_Pi','; q2^{2} (MeV^{2}); E_{l} (MeV^{2}); Mmiss^{2} (MeV^{2})',4,-2E6,14E6,10,0,2600,10,-2E6,14E6)

for i in range(t.GetEntries()):
    t.GetEntry(i)
    h_El.Fill(t.FitVar_El_mLc,t.sw_sig)
    h_q2.Fill(t.FitVar_q2_mLc,t.sw_sig)
    h_Mmiss.Fill(t.FitVar_Mmiss2_mLc,t.sw_sig)
    h_mLc.Fill(t.Lc_M,t.sw_sig)
    h_dataSS.Fill(t.FitVar_q2_mLc, t.FitVar_El_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig)

for i in range(t1.GetEntries()):
    t1.GetEntry(i)
    h_El.Fill(t1.FitVar_El_mLc,t1.sw_sig)
    h_q2.Fill(t1.FitVar_q2_mLc,t1.sw_sig)
    h_Mmiss.Fill(t1.FitVar_Mmiss2_mLc,t1.sw_sig)
    h_mLc.Fill(t1.Lc_M,t1.sw_sig)
    h_dataSS.Fill(t1.FitVar_q2_mLc, t1.FitVar_El_mLc, t1.FitVar_Mmiss2_mLc, t1.sw_sig)
    

f_Pi_magUp = r.TFile('/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MISID/SameSign/Pi_sample_MagUp_Kenr_sw_withCF.root','READ')
t2 = f_Pi_magUp.Get('DecayTree')
f_Pi_magDown = r.TFile('/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MISID/SameSign/Pi_sample_MagDown_Kenr_sw_withCF.root','READ')
t3 = f_Pi_magDown.Get('DecayTree')

for i in range(t2.GetEntries()):
    t2.GetEntry(i)
    h1_El.Fill(t2.FitVar_El_mLc,t2.sw_sig*t2.w_recomu_CF*10)
    h1_q2.Fill(t2.FitVar_q2_mLc,t2.sw_sig*t2.w_recomu_CF*10)
    h1_Mmiss.Fill(t2.FitVar_Mmiss2_mLc,t2.sw_sig*t2.w_recomu_CF*10)
    h1_mLc.Fill(t2.Lc_M,t2.sw_sig*t2.w_recomu_CF*10)
    h_Pi.Fill(t2.FitVar_q2_mLc, t2.FitVar_El_mLc, t2.FitVar_Mmiss2_mLc, t2.sw_sig*t2.w_recomu_CF*10)

for i in range(t3.GetEntries()):
    t3.GetEntry(i)
    h1_El.Fill(t3.FitVar_El_mLc,t3.sw_sig*t3.w_recomu_CF*10)
    h1_q2.Fill(t3.FitVar_q2_mLc,t3.sw_sig*t3.w_recomu_CF*10)
    h1_Mmiss.Fill(t3.FitVar_Mmiss2_mLc,t3.sw_sig*t3.w_recomu_CF*10)
    h1_mLc.Fill(t3.Lc_M,t3.sw_sig*t3.w_recomu_CF*10)
    h_Pi.Fill(t3.FitVar_q2_mLc, t3.FitVar_El_mLc, t3.FitVar_Mmiss2_mLc, t3.sw_sig*t3.w_recomu_CF*10)

print(h_mLc.Integral(), h1_mLc.Integral())

for i in range(h_dataSS.GetNbinsX()):
    for j in range(h_dataSS.GetNbinsY()):
        for k in range(h_dataSS.GetNbinsZ()):
            ndata = h_dataSS.GetBinContent(i,j,k)
            nPi =  h_Pi.GetBinContent(i,j,k)
            '''
            if ndata<1:
                ndata=0.
                print('Ndata rounded')
            if nPi<1:
                nPi=0
                print('NPi rounded')
            '''
            if ndata!=0:
                fPi = nPi*1./ndata
            else:
                fPi=0.
            print(ndata, nPi, fPi)

c = r.TCanvas('c','',1500,500)
c.Divide(3,1)
c.cd(1)
h_El.Draw()
h1_El.SetLineColor(r.kBlue)
h1_El.Draw('same')
c.cd(2)
h_q2.Draw()
h1_q2.SetLineColor(r.kBlue)
h1_q2.Draw('same')
c.cd(3)
h_Mmiss.Draw()
h1_Mmiss.SetLineColor(r.kBlue)
h1_Mmiss.Draw('same')

'''
c1 = r.TCanvas('c1','',1000,1000)
c1.Divide(2,2)
c1.cd(1)
'''


