import ROOT as r

datadir = '/disk/lhcb_data2/RLcMuonic2016/Data/'
polarities = ['MagUp','MagDown']

h_new = r.TH3F('h_new','; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
h_old = r.TH3F('h_old','; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)

for polarity in polarities:
    fnameNew = datadir+'Lb_Data_'+polarity+'_preselected_Kenr_sw.root'
    fnew = r.TFile(fnameNew,'READ')
    tnew = fnew.Get('DecayTree')

    for i in range(tnew.GetEntries()):
        tnew.GetEntry(i)
        h_new.Fill(tnew.FitVar_q2_mLc/1E6,tnew.FitVar_El_mLc, tnew.FitVar_Mmiss2_mLc/1E6,tnew.sw_sig)
    
    fnameOld = datadir+'/KenrWithGhosts/Lb_Data_'+polarity+'_preselected_Kenr_sw.root'
    fold = r.TFile(fnameOld,'READ')
    told = fold.Get('DecayTree')

    for i in range(told.GetEntries()):
        told.GetEntry(i)
        h_old.Fill(told.FitVar_q2_mLc/1E6,told.FitVar_El_mLc, told.FitVar_Mmiss2_mLc/1E6,told.sw_sig)


print('Data Yield new Kenriched selection: {:.2f}'.format(h_new.Integral()))
print('Data Yield old Kenriched selection: {:.2f}'.format(h_old.Integral()))
h_q2_new = h_new.ProjectionX()
h_q2_old = h_old.ProjectionX()
h_El_new = h_new.ProjectionY()
h_El_old = h_old.ProjectionY()
h_Mmiss_new = h_new.ProjectionZ()
h_Mmiss_old = h_old.ProjectionZ()

c = r.TCanvas('c','Templates',1500,500)
c.Divide(3,1)
c.cd(1)
h_q2_old.SetMarkerStyle(20)
h_q2_old.SetMarkerSize(1)
h_q2_old.SetLineColor(r.kRed+1)
h_q2_old.SetMarkerColor(r.kRed+1)
h_q2_old.Draw('PE0')
h_q2_new.SetMarkerStyle(20)
h_q2_new.SetMarkerSize(1)
h_q2_new.SetLineColor(r.kBlack)
h_q2_new.Draw('PE0 sames')
l = r.TLegend(0.1,0.75,0.4,0.9)
l.AddEntry(h_q2_new,'New selection','pl')
l.AddEntry(h_q2_old,'Old selection','pl')
l.Draw()
c.cd(2)
h_El_old.SetMarkerStyle(20)
h_El_old.SetMarkerSize(1)
h_El_old.SetLineColor(r.kRed+1)
h_El_old.SetMarkerColor(r.kRed+1)
h_El_old.Draw('PE0')
h_El_new.SetMarkerStyle(20)
h_El_new.SetMarkerSize(1)
h_El_new.SetLineColor(r.kBlack)
h_El_new.Draw('PE0 sames')
c.cd(3)
h_Mmiss_old.SetMarkerStyle(20)
h_Mmiss_old.SetMarkerSize(1)
h_Mmiss_old.SetLineColor(r.kRed+1)
h_Mmiss_old.SetMarkerColor(r.kRed+1)
h_Mmiss_old.Draw('PE0')
h_Mmiss_new.SetMarkerStyle(20)
h_Mmiss_new.SetMarkerSize(1)
h_Mmiss_new.SetLineColor(r.kBlack)
h_Mmiss_new.Draw('PE0 sames')

