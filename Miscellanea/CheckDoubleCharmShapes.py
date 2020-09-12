import ROOT as r

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
polarities = ['MagUp','MagDown']
categories = ['iso','Kenr']
h_LcDs_iso = r.TH3F('h_LcDs_iso','h_LcDs_iso;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)
h_Lc2625Ds_iso = r.TH3F('h_Lc2625Ds_iso','h_Lc2625Ds_iso;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)
h_Lc2593Ds_iso = r.TH3F('h_Lc2593Ds_iso','h_Lc2593Ds_iso;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)
h_LcDs_Kenr = r.TH3F('h_LcDs_Kenr','h_LcDs_Kenr;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)
h_Lc2625Ds_Kenr = r.TH3F('h_Lc2625Ds_Kenr','h_Lc2625Ds_Kenr;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)
h_Lc2593Ds_Kenr = r.TH3F('h_Lc2593Ds_Kenr','h_Lc2593Ds_Kenr;q^{2};E_{l};M_{miss}^2',4,-2,14,10,0,2600,10,-2,14)

h_LcDs = {'iso':h_LcDs_iso,'Kenr':h_LcDs_Kenr}
h_Lc2625Ds = {'iso':h_Lc2625Ds_iso,'Kenr':h_Lc2625Ds_Kenr}
h_Lc2593Ds = {'iso':h_Lc2593Ds_iso,'Kenr':h_Lc2593Ds_Kenr}

for category in categories:
    for polarity in polarities:
        fname = 'Lb_LcDs_'+polarity+'_full_preselected_'+category+'.root'
        f = r.TFile(filedir+fname,'Read')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            weight = t.Event_PIDCalibEffWeight
            h_LcDs[category].Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)

        fname = 'Lb_Lc2625Ds_'+polarity+'_full_preselected_'+category+'.root'
        f = r.TFile(filedir+fname,'Read')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            weight = t.Event_PIDCalibEffWeight
            h_Lc2625Ds[category].Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
        
        fname = 'Lb_Lc2593Ds_'+polarity+'_full_preselected_'+category+'.root'
        f = r.TFile(filedir+fname,'Read')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            weight = t.Event_PIDCalibEffWeight
            h_Lc2593Ds[category].Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)

h_LcDs_iso_q2 = h_LcDs_iso.ProjectionX()
h_LcDs_iso_El = h_LcDs_iso.ProjectionY()
h_LcDs_iso_Mmiss = h_LcDs_iso.ProjectionZ()

h_Lc2625Ds_iso_q2 = h_Lc2625Ds_iso.ProjectionX()
h_Lc2625Ds_iso_El = h_Lc2625Ds_iso.ProjectionY()
h_Lc2625Ds_iso_Mmiss = h_Lc2625Ds_iso.ProjectionZ()

h_Lc2593Ds_iso_q2 = h_Lc2593Ds_iso.ProjectionX()
h_Lc2593Ds_iso_El = h_Lc2593Ds_iso.ProjectionY()
h_Lc2593Ds_iso_Mmiss = h_Lc2593Ds_iso.ProjectionZ()

h_LcDs_Kenr_q2 = h_LcDs_Kenr.ProjectionX()
h_LcDs_Kenr_El = h_LcDs_Kenr.ProjectionY()
h_LcDs_Kenr_Mmiss = h_LcDs_Kenr.ProjectionZ()

h_Lc2625Ds_Kenr_q2 = h_Lc2625Ds_Kenr.ProjectionX()
h_Lc2625Ds_Kenr_El = h_Lc2625Ds_Kenr.ProjectionY()
h_Lc2625Ds_Kenr_Mmiss = h_Lc2625Ds_Kenr.ProjectionZ()

h_Lc2593Ds_Kenr_q2 = h_Lc2593Ds_Kenr.ProjectionX()
h_Lc2593Ds_Kenr_El = h_Lc2593Ds_Kenr.ProjectionY()
h_Lc2593Ds_Kenr_Mmiss = h_Lc2593Ds_Kenr.ProjectionZ()

'''
h_LcDs_iso_q2 = ScaleHisto(h_LcDs_iso_q2,1)
h_LcDs_iso_El = ScaleHisto(h_LcDs_iso_El,1)
h_LcDs_iso_Mmiss = ScaleHisto(h_LcDs_iso_Mmiss,1)

h_LcDs_Kenr_q2 = ScaleHisto(h_LcDs_Kenr_q2,1)
h_LcDs_Kenr_El = ScaleHisto(h_LcDs_Kenr_El,1)
h_LcDs_Kenr_Mmiss = ScaleHisto(h_LcDs_Kenr_Mmiss,1)

h_Lc2625Ds_iso_q2 = ScaleHisto(h_Lc2625Ds_iso_q2,1)
h_Lc2625Ds_iso_El = ScaleHisto(h_Lc2625Ds_iso_El,1)
h_Lc2625Ds_iso_Mmiss = ScaleHisto(h_Lc2625Ds_iso_Mmiss,1)

h_Lc2625Ds_Kenr_q2 = ScaleHisto(h_Lc2625Ds_Kenr_q2,1)
h_Lc2625Ds_Kenr_El = ScaleHisto(h_Lc2625Ds_Kenr_El,1)
h_Lc2625Ds_Kenr_Mmiss = ScaleHisto(h_Lc2625Ds_Kenr_Mmiss,1)

h_Lc2593Ds_iso_q2 = ScaleHisto(h_Lc2593Ds_iso_q2,1)
h_Lc2593Ds_iso_El = ScaleHisto(h_Lc2593Ds_iso_El,1)
h_Lc2593Ds_iso_Mmiss = ScaleHisto(h_Lc2593Ds_iso_Mmiss,1)

h_Lc2593Ds_Kenr_q2 = ScaleHisto(h_Lc2593Ds_Kenr_q2,1)
h_Lc2593Ds_Kenr_El = ScaleHisto(h_Lc2593Ds_Kenr_El,1)
h_Lc2593Ds_Kenr_Mmiss = ScaleHisto(h_Lc2593Ds_Kenr_Mmiss,1)
'''

c = r.TCanvas('c','c',1500,500)
c.Divide(3,1)
c.cd(1)
#h_LcDs_iso_q2.GetYaxis().SetRangeUser(0,0.7)
h_LcDs_iso_q2.Draw('hist')
h_Lc2625Ds_iso_q2.SetLineColor(r.kRed)
h_Lc2625Ds_iso_q2.Draw('hist sames')
h_Lc2593Ds_iso_q2.SetLineColor(r.kBlack)
h_Lc2593Ds_iso_q2.Draw('hist sames')
l = r.TLegend(0.1,0.75,0.4,0.9)
l.AddEntry(h_LcDs_iso_q2,'#Lambda_{c}X_{c}','pl')
l.AddEntry(h_Lc2625Ds_iso_q2,'#Lambda_{c}^{*}(2625)X_{c}','pl')
l.AddEntry(h_Lc2593Ds_iso_q2,'#Lambda_{c}^{*}(2593)X_{c}','pl')
l.Draw()
c.cd(2)
#h_LcDs_iso_El.GetYaxis().SetRangeUser(0,0.7)
h_LcDs_iso_El.Draw('hist')
h_Lc2625Ds_iso_El.SetLineColor(r.kRed+1)
h_Lc2625Ds_iso_El.Draw('hist sames')
h_Lc2593Ds_iso_El.SetLineColor(r.kBlack)
h_Lc2593Ds_iso_El.Draw('hist sames')
c.cd(3)
#h_LcDs_iso_Mmiss.GetYaxis().SetRangeUser(0,0.7)
h_LcDs_iso_Mmiss.Draw('hist')
h_Lc2625Ds_iso_Mmiss.SetLineColor(r.kRed+1)
h_Lc2625Ds_iso_Mmiss.Draw('hist sames')
h_Lc2593Ds_iso_Mmiss.SetLineColor(r.kBlack)
h_Lc2593Ds_iso_Mmiss.Draw('hist sames')

c1 = r.TCanvas('c1','c1',1500,500)
c1.Divide(3,1)
c1.cd(1)
#h_LcDs_Kenr_q2.GetYaxis().SetRangeUser(0,1)
h_LcDs_Kenr_q2.Draw('hist')
h_Lc2625Ds_Kenr_q2.SetLineColor(r.kRed+1)
h_Lc2625Ds_Kenr_q2.Draw('hist sames')
h_Lc2593Ds_Kenr_q2.SetLineColor(r.kBlack)
h_Lc2593Ds_Kenr_q2.Draw('hist sames')
l1 = r.TLegend(0.1,0.75,0.4,0.9)
l1.AddEntry(h_LcDs_iso_q2,'#Lambda_{c}X_{c}','pl')
l1.AddEntry(h_Lc2625Ds_iso_q2,'#Lambda_{c}^{*}(2625)X_{c}','pl')
l1.AddEntry(h_Lc2593Ds_iso_q2,'#Lambda_{c}^{*}(2593)X_{c}','pl')
l1.Draw()
c1.cd(2)
#h_LcDs_Kenr_El.GetYaxis().SetRangeUser(0,1)
h_LcDs_Kenr_El.Draw('hist')
h_Lc2625Ds_Kenr_El.SetLineColor(r.kRed+1)
h_Lc2625Ds_Kenr_El.Draw('hist sames')
h_Lc2593Ds_Kenr_El.SetLineColor(r.kBlack)
h_Lc2593Ds_Kenr_El.Draw('hist sames')
c1.cd(3)
#h_LcDs_Kenr_Mmiss.GetYaxis().SetRangeUser(0,1)
h_LcDs_Kenr_Mmiss.Draw('hist')
h_Lc2625Ds_Kenr_Mmiss.SetLineColor(r.kRed+1)
h_Lc2625Ds_Kenr_Mmiss.Draw('hist sames')
h_Lc2593Ds_Kenr_Mmiss.SetLineColor(r.kBlack)
h_Lc2593Ds_Kenr_Mmiss.Draw('hist sames')
