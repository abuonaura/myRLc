'''
Author: Annarita Buonaura
Data: 5 April 2019
Purpose: Draw fit templates of dataSS and misidSS samples for the different isolation cases (no iso, iso applied, anti-isolated sample)
+ Draw fit templates of data and (dataSS) combinatorial bkg
'''

import ROOT as r
import os, sys, getopt

datadir = '/disk/lhcb_data2/RLcMuonic2016/'
polarities=['MagUp','MagDown']
#polarities=['MagDown']
particles=['K','Pi']
sample_suffix = {'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root','Lcpipi':'_Lcpipi_sw.root'}
sample_suffix_misid = {'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root','Lcpipi':'_Lcpipi_sw_withCF.root'}
sample_suffix_comb = {'iso':'_sw_noMISID_iso.root','Kenriched':'_sw_noMISID_Kenr.root','Lcpipi':'_sw_noMISID_Lcpipi.root'}

def makeTemplateDataSS(polarity,sample):
    histo = r.TH3F('h_'+polarity+'SS','; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    datafname = datadir + 'Data/Lb_DataSS_'+polarity+'_preselected'+sample_suffix[sample]

    dataf = r.TFile(datafname, 'Read')
    datat = dataf.Get('DecayTree')
    print(datat.GetEntries())
    for i in range(datat.GetEntries()):
        datat.GetEntry(i)
        histo.Fill(datat.FitVar_q2_mLc/1E6,datat.FitVar_El_mLc, datat.FitVar_Mmiss2_mLc/1E6,datat.sw_sig)
    print(histo.Integral(), histo.GetEntries())
    return histo

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

def makeTemplateCombinatorial(polarity,sample):
    histo = r.TH3F('h_'+polarity+'Comb','; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    datafname = datadir + 'Data/Lb_DataSS_'+polarity+sample_suffix_comb[sample]
    dataf = r.TFile(datafname, 'Read')
    datat = dataf.Get('DecayTree')
    print(datat.GetEntries())
    for i in range(datat.GetEntries()):
        datat.GetEntry(i)
        if datat.Lb_M<5620:
            histo.Fill(datat.FitVar_q2_mLc/1E6,datat.FitVar_El_mLc, datat.FitVar_Mmiss2_mLc/1E6,datat.sw_sig*datat.w_MISID)
    print(histo.Integral(), histo.GetEntries())
    return histo


def makeTemplateMISIDSS(particle,polarity,sample):
    histo = r.TH3F('h_'+particle+'_'+polarity,'; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    fname = datadir + 'MISID/SameSign/'+particle+'_sample_'+polarity+sample_suffix_misid[sample]
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
    opts, args = getopt.getopt(sys.argv[1:], "",["Lcpipi","iso","Kenriched","norm"])
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
    
    #Plotting the DataSS templates and the MISID SS templates
    h_MagUpSS = makeTemplateDataSS('MagUp',sample)
    h_MagDownSS = makeTemplateDataSS('MagDown',sample)
    h_dataSS = r.TH3F('h_dataSS','3D histo Data; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_dataSS.Add(h_MagUpSS, h_MagDownSS)
    
    h_Pi = r.TH3F('h_Pi','3D histo MISID pi; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_PiUp = makeTemplateMISIDSS('Pi','MagUp',sample)
    h_PiDown = makeTemplateMISIDSS('Pi','MagDown',sample)
    h_Pi.Add(h_PiUp, h_PiDown)
    
    h_K = r.TH3F('h_K','3D histo MISID K; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_KUp = makeTemplateMISIDSS('K','MagUp',sample)
    h_KDown = makeTemplateMISIDSS('K','MagDown',sample)
    h_K.Add(h_KUp, h_KDown)


    h_q2_dataSS = h_dataSS.ProjectionX()
    h_El_dataSS = h_dataSS.ProjectionY()
    h_Mmiss_dataSS = h_dataSS.ProjectionZ()
    h_q2_Pi = h_Pi.ProjectionX()
    h_El_Pi = h_Pi.ProjectionY()
    h_Mmiss_Pi = h_Pi.ProjectionZ()
    h_q2_K = h_K.ProjectionX()
    h_El_K = h_K.ProjectionY()
    h_Mmiss_K = h_K.ProjectionZ()
    if normalise:
        h_q2_dataSS = ScaleHisto(h_q2_dataSS,1)
        h_El_dataSS = ScaleHisto(h_El_dataSS,1)
        h_Mmiss_dataSS = ScaleHisto(h_Mmiss_dataSS,1)
        h_q2_Pi = ScaleHisto(h_q2_Pi,1)
        h_El_Pi = ScaleHisto(h_El_Pi,1)
        h_Mmiss_Pi = ScaleHisto(h_Mmiss_Pi,1)
        h_q2_K = ScaleHisto(h_q2_K,1)
        h_El_K = ScaleHisto(h_El_K,1)
        h_Mmiss_K = ScaleHisto(h_Mmiss_K,1)

    r.gStyle.SetOptStat(0)
    c = r.TCanvas('c','Templates',1500,500)
    c.Divide(3,1)
    c.cd(1)
    if normalise:
        h_q2_dataSS.GetYaxis().SetRangeUser(0,0.5)
    h_q2_dataSS.SetMarkerStyle(20)
    h_q2_dataSS.SetMarkerSize(1)
    h_q2_dataSS.SetLineColor(r.kBlack)
    h_q2_dataSS.Draw('PE0')
    h_q2_Pi.SetLineColor(r.kBlue)
    h_q2_Pi.Draw('histo sames')
    h_q2_K.SetLineColor(r.kRed)
    h_q2_K.Draw('histo sames')
    l = r.TLegend(0.1,0.75,0.4,0.9)
    l.AddEntry(h_q2_dataSS,'Data SS','pl')
    l.AddEntry(h_q2_Pi,'MISID Pi Data SS','l')
    l.AddEntry(h_q2_K,'MISID K Data SS','l')
    l.Draw()
    c.cd(2)
    if normalise:
        h_El_dataSS.GetYaxis().SetRangeUser(0,0.3)
    h_El_dataSS.SetMarkerStyle(20)
    h_El_dataSS.SetMarkerSize(1)
    h_El_dataSS.SetLineColor(r.kBlack)
    h_El_dataSS.Draw('PE0')
    h_El_Pi.SetLineColor(r.kBlue)
    h_El_Pi.Draw('histo sames')
    h_El_K.SetLineColor(r.kRed)
    h_El_K.Draw('histo sames')
    c.cd(3)
    if normalise:
        h_Mmiss_dataSS.GetYaxis().SetRangeUser(0,0.3)
    h_Mmiss_dataSS.SetMarkerStyle(20)
    h_Mmiss_dataSS.SetMarkerSize(1)
    h_Mmiss_dataSS.SetLineColor(r.kBlack)
    h_Mmiss_dataSS.Draw('PE0')
    h_Mmiss_Pi.SetLineColor(r.kBlue)
    h_Mmiss_Pi.Draw('histo sames')
    h_Mmiss_K.SetLineColor(r.kRed)
    h_Mmiss_K.Draw('histo sames')
    if normalise:
        c.SaveAs('plots/TemplatesSS_'+sample+'_normalised.C')
    else:
        c.SaveAs('plots/TemplatesSS_'+sample+'.C')

    print('DataSS Yield: {:.2f}'.format(h_dataSS.Integral()))
    print('DataSS pi MISID Yield: {:.2f}'.format(h_Pi.Integral()))
    print('    fraction: {:.2f}'.format(h_Pi.Integral()*100/h_dataSS.Integral()))
    print('DataSS K MISID Yield: {:.2f}'.format(h_K.Integral()))
    print('    fraction: {:.2f}'.format(h_K.Integral()*100/h_dataSS.Integral()))

    #Plotting data and combinatorial bkg templates
    h_MagUp = makeTemplateData('MagUp',sample)
    h_MagDown = makeTemplateData('MagDown',sample)
    h_data = r.TH3F('h_data','3D histo Data; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_data.Add(h_MagUp, h_MagDown)
    
    h_MagUpComb = makeTemplateCombinatorial('MagUp',sample)
    h_MagDownComb = makeTemplateCombinatorial('MagDown',sample)
    h_comb = r.TH3F('h_comb','3D histo Comb bkg; q^{2} (Gev^{2}); E_{l} (MeV^{2}); M_{miss}^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_comb.Add(h_MagUpComb, h_MagDownComb)

    h_q2_data = h_data.ProjectionX()
    h_El_data = h_data.ProjectionY()
    h_Mmiss_data = h_data.ProjectionZ()
    h_q2_comb = h_comb.ProjectionX()
    h_El_comb = h_comb.ProjectionY()
    h_Mmiss_comb = h_comb.ProjectionZ()
    if normalise:
        h_q2_data = ScaleHisto(h_q2_data,1)
        h_El_data = ScaleHisto(h_El_data,1)
        h_Mmiss_data = ScaleHisto(h_Mmiss_data,1)
        h_q2_comb = ScaleHisto(h_q2_comb,1)
        h_El_comb = ScaleHisto(h_El_comb,1)
        h_Mmiss_comb = ScaleHisto(h_Mmiss_comb,1)

    c1 = r.TCanvas('c1','Combinatorial BKg',1500,500)
    c1.Divide(3,1)
    c1.cd(1)
    h_q2_data.SetMarkerStyle(20)
    h_q2_data.SetMarkerSize(1)
    h_q2_data.SetLineColor(r.kBlack)
    h_q2_data.GetYaxis().SetLabelSize(0.045)
    h_q2_data.GetXaxis().SetLabelSize(0.045)
    h_q2_data.GetXaxis().SetTitleSize(0.045)
    h_q2_data.Draw('PE0')
    h_q2_comb.SetLineColor(r.kBlue)
    h_q2_comb.Draw('histo sames')
    l1 = r.TLegend(0.1,0.75,0.4,0.9)
    l1.AddEntry(h_q2_data,'Data ','pl')
    l1.AddEntry(h_q2_comb,'CombBkg','l')
    l1.Draw()
    c1.cd(2)
    h_El_data.SetMarkerStyle(20)
    h_El_data.SetMarkerSize(1)
    h_El_data.SetLineColor(r.kBlack)
    h_El_data.GetYaxis().SetLabelSize(0.045)
    h_El_data.GetXaxis().SetLabelSize(0.045)
    h_El_data.GetXaxis().SetTitleSize(0.045)
    h_El_data.Draw('PE0')
    h_El_comb.SetLineColor(r.kBlue)
    h_El_comb.Draw('histo sames')
    c1.cd(3)
    h_Mmiss_data.SetMarkerStyle(20)
    h_Mmiss_data.SetMarkerSize(1)
    h_Mmiss_data.SetLineColor(r.kBlack)
    h_Mmiss_data.GetYaxis().SetLabelSize(0.045)
    h_Mmiss_data.GetXaxis().SetLabelSize(0.045)
    h_Mmiss_data.GetXaxis().SetTitleSize(0.045)
    h_Mmiss_data.Draw('PE0')
    h_Mmiss_comb.SetLineColor(r.kBlue)
    h_Mmiss_comb.Draw('histo sames')
    if normalise:
        c1.SaveAs('plots/TemplatesComb_'+sample+'_normalised.C')
    else:
        c1.SaveAs('plots/TemplatesComb_'+sample+'.C')
    
    print('Data Yield: {:.2f}'.format(h_data.Integral()))
    print('Comb Bkg Yield: {:.2f}'.format(h_comb.Integral()))
    print('    fraction: {:.2f}'.format(h_comb.Integral()*100/h_data.Integral()))

