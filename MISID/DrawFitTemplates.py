'''
Author: Annarita Buonaura
Data: 5 April 2019
Purpose: Draw fit templates of data and misid samples for the different isolation cases (no iso, iso applied, anti-isolated sample)
'''

import ROOT as r
import os, sys, getopt

datadir = '/disk/lhcb_data2/RLcMuonic2016/'

polarities=['MagUp','MagDown']
#polarities=['MagDown']
particles=['K','Pi']
sample_suffix = {'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root','Lcpipi':'_Lcpipi_sw.root'}
sample_suffix_misid = {'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root','Lcpipi':'_Lcpipi_sw_withCF.root'}

def makeTemplateData(polarity,sample):
    histo = r.TH3F('h_'+polarity,'; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    datafname = datadir + 'Data/Lb_Data_'+polarity+'_preselected'+sample_suffix[sample]
    dataf = r.TFile(datafname, 'Read')
    datat = dataf.Get('DecayTree')
    print(datat.GetEntries())
    for i in range(datat.GetEntries()):
        datat.GetEntry(i)
        histo.Fill(datat.FitVar_q2_mLc/1E6,datat.FitVar_El_mLc, datat.FitVar_Mmiss2_mLc/1E6,datat.sw_sig)
    print(histo.Integral(), histo.GetEntries())
    return histo

def makeTemplateMISID(particle,polarity,sample):
    histo = r.TH3F('h_'+particle+'_'+polarity,'; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    fname = datadir + 'MISID/OppositeSign/'+particle+'_sample_'+polarity+sample_suffix_misid[sample]
    f = r.TFile(fname, 'Read')
    t = f.Get('DecayTree')
    print(t.GetEntries())
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        histo.Fill(t.FitVar_q2_mLc/1E6,t.FitVar_El_mLc, t.FitVar_Mmiss2_mLc/1E6,t.sw_sig*t.w_recomu_CF*10)
    print(histo.Integral(), histo.GetEntries())
    return histo

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h



if __name__ == '__main__':

    normalise = False
    opts, args = getopt.getopt(sys.argv[1:], "",["iso","Kenriched","Lcpipi","norm"])
    print (opts,args)
    for o, a in opts:
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'
        if o in ("--Lcpipi",):
            sample = 'Lcpipi'
        if o in ("--norm",):
            normalise = True
    
    h_MagUp = makeTemplateData('MagUp',sample)
    h_MagDown = makeTemplateData('MagDown',sample)
    h_data = r.TH3F('h_data','3D histo Data; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_data.Add(h_MagUp, h_MagDown)
    
    h_Pi = r.TH3F('h_Pi','3D histo MISID pi; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_PiUp = makeTemplateMISID('Pi','MagUp',sample)
    h_PiDown = makeTemplateMISID('Pi','MagDown',sample)
    h_Pi.Add(h_PiUp, h_PiDown)
    
    h_K = r.TH3F('h_K','3D histo MISID K; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mi_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_KUp = makeTemplateMISID('K','MagUp',sample)
    h_KDown = makeTemplateMISID('K','MagDown',sample)
    h_K.Add(h_KUp, h_KDown)

    h_q2_data = h_data.ProjectionX()
    h_q2_Pi = h_Pi.ProjectionX()
    h_q2_K = h_K.ProjectionX()
    h_El_data = h_data.ProjectionY()
    h_El_Pi = h_Pi.ProjectionY()
    h_El_K = h_K.ProjectionY()
    h_Mmiss_Pi = h_Pi.ProjectionZ()
    h_Mmiss_data = h_data.ProjectionZ()
    h_Mmiss_K = h_K.ProjectionZ()
    if normalise:
        h_q2_data = ScaleHisto(h_q2_data,1)
        h_El_data = ScaleHisto(h_El_data,1)
        h_Mmiss_data = ScaleHisto(h_Mmiss_data,1)
        h_q2_K = ScaleHisto(h_q2_K,1)
        h_El_K = ScaleHisto(h_El_K,1)
        h_Mmiss_K = ScaleHisto(h_Mmiss_K,1)
        h_q2_Pi = ScaleHisto(h_q2_Pi,1)
        h_El_Pi = ScaleHisto(h_El_Pi,1)
        h_Mmiss_Pi = ScaleHisto(h_Mmiss_Pi,1)

    r.gStyle.SetOptStat(0)
    c = r.TCanvas('c','Templates',1500,500)
    c.Divide(3,1)
    c.cd(1)
    h_q2_data.SetMarkerStyle(20)
    h_q2_data.SetMarkerSize(1)
    h_q2_data.SetLineColor(r.kBlack)
    h_q2_data.GetYaxis().SetLabelSize(0.045)
    h_q2_data.GetXaxis().SetLabelSize(0.045)
    h_q2_data.GetXaxis().SetTitleSize(0.045)
    h_q2_data.Draw('PE0')
    h_q2_Pi.SetLineColor(r.kBlue)
    h_q2_Pi.Draw('histo sames')
    h_q2_K.SetLineColor(r.kRed)
    h_q2_K.Draw('histo sames')
    l = r.TLegend(0.1,0.75,0.4,0.9)
    l.AddEntry(h_q2_data,'Data','pl')
    l.AddEntry(h_q2_Pi,'MISID #pi ','l')
    l.AddEntry(h_q2_K,'MISID K','l')
    l.Draw()
    c.cd(2)
    h_El_data.SetMarkerStyle(20)
    h_El_data.SetMarkerSize(1)
    h_El_data.SetLineColor(r.kBlack)
    if normalise:
        h_El_data.GetYaxis().SetRangeUser(0,0.3)
    h_El_data.GetYaxis().SetLabelSize(0.045)
    h_El_data.GetXaxis().SetLabelSize(0.045)
    h_El_data.GetXaxis().SetTitleSize(0.045)
    h_El_data.Draw('PE0')
    h_El_Pi.SetLineColor(r.kBlue)
    h_El_Pi.Draw('histo sames')
    h_El_K.SetLineColor(r.kRed)
    h_El_K.Draw('histo sames')
    c.cd(3)
    h_Mmiss_data.SetMarkerStyle(20)
    h_Mmiss_data.SetMarkerSize(1)
    h_Mmiss_data.SetLineColor(r.kBlack)
    h_Mmiss_data.GetYaxis().SetLabelSize(0.045)
    h_Mmiss_data.GetXaxis().SetLabelSize(0.045)
    h_Mmiss_data.GetXaxis().SetTitleSize(0.045)
    h_Mmiss_data.Draw('PE0')
    h_Mmiss_Pi.SetLineColor(r.kBlue)
    h_Mmiss_Pi.Draw('histo sames')
    h_Mmiss_K.SetLineColor(r.kRed)
    h_Mmiss_K.Draw('histo sames')
    if normalise:
        c.SaveAs('plots/Templates_'+sample+'_normalised.C')
    else:
        c.SaveAs('plots/Templates_'+sample+'.C')

    print('Data Yield: {:.2f}'.format(h_data.Integral()))
    print('Data pi MISID Yield: {:.2f}'.format(h_Pi.Integral()))
    print('    fraction: {:.2f}'.format(h_Pi.Integral()*100/h_data.Integral()))
    print('Data K MISID Yield: {:.2f}'.format(h_K.Integral()))
    print('    fraction: {:.2f}'.format(h_K.Integral()*100/h_data.Integral()))
'''
sample = 'antiso'
h = r.TH3F('h','3D histo Data; Mmiss^{2} (Gev^{2}); E_{l} (MeV^{2}); q^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
for polarity in polarities:
    datafname = datadir + 'Data/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
    dataf = r.TFile(datafname, 'Read')
    datat = dataf.Get('DecayTree')
    for i in range(datat.GetEntries()):
        datat.GetEntry(i)
        h.Fill(datat.FitVar_q2_mLc/1E6,datat.FitVar_El_mLc, datat.FitVar_Mmiss2_mLc/1E6,datat.sw_sig)
c = r.TCanvas()
h.Draw('box')
'''
