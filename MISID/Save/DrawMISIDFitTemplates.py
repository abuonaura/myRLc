'''
Author: Annarita Buonaura
Data: 5 April 2019
Purpose: Draw fit templates of data and misid samples for the different isolation cases (no iso, iso applied, anti-isolated sample)
'''

import ROOT as r
import os, sys, getopt

datadir = '$FILEDIR/'

#polarities=['MagUp','MagDown']
polarities=['MagDown']
particles=['K','Pi']
sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
sample_suffix_misid = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}


def makeTemplateMISID(particle,polarity,sample):
    histo = r.TH3F('h_'+particle+'_'+polarity,'; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    fname = datadir + 'MISID/OppositeSign/'+particle+'_sample_'+polarity+sample_suffix_misid[sample]
    f = r.TFile(fname, 'Read')
    t = f.Get('DecayTree')
    print(t.GetEntries())
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        histo.Fill(t.FitVar_q2_mLc/1E6,t.FitVar_El_mLc, t.FitVar_Mmiss2_mLc/1E6,t.sw_sig*t.w_recomu_CF*10)
    print(histo.Integral(), histo.GetEntries())
    return histo


if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], "",["full","iso","Kenriched"])
    print (opts,args)
    for o, a in opts:
        if o in ("--full",):
            sample = 'full'
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'
    
    
    h_Pi = r.TH3F('h_Pi','3D histo MISID pi; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_PiUp = makeTemplateMISID('Pi','MagUp',sample)
    h_PiDown = makeTemplateMISID('Pi','MagDown',sample)
    h_Pi.Add(h_PiUp, h_PiDown)
    
    h_K = r.TH3F('h_K','3D histo MISID K; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_KUp = makeTemplateMISID('K','MagUp',sample)
    h_KDown = makeTemplateMISID('K','MagDown',sample)
    h_K.Add(h_KUp, h_KDown)

    r.gStyle.SetOptStat(0)
    c = r.TCanvas('c','Templates',1600,500)
    c.Divide(3,1)
    c.cd(1)
    h_q2_Pi = h_Pi.ProjectionX()
    h_q2_Pi.SetTitle("")
    h_q2_Pi.GetYaxis().SetRangeUser(0,1.6E4)
    h_q2_Pi.GetYaxis().SetLabelSize(0.045)
    h_q2_Pi.GetXaxis().SetLabelSize(0.045)
    h_q2_Pi.GetXaxis().SetTitleSize(0.045)
    h_q2_Pi.SetLineColor(r.kBlue)
    h_q2_Pi.Draw('histo')
    h_q2_K = h_K.ProjectionX()
    h_q2_K.SetLineColor(r.kRed)
    h_q2_K.Draw('histo sames')
    l = r.TLegend(0.1,0.75,0.4,0.9)
    l.AddEntry(h_q2_Pi,'#pi sample','l')
    l.AddEntry(h_q2_K,'K sample','l')
    l.Draw()
    c.cd(2)
    h_El_Pi = h_Pi.ProjectionY()
    h_El_Pi.SetTitle("")
    h_El_Pi.GetYaxis().SetRangeUser(0,7.5E3)
    h_El_Pi.GetYaxis().SetLabelSize(0.045)
    h_El_Pi.GetXaxis().SetLabelSize(0.045)
    h_El_Pi.GetXaxis().SetTitleSize(0.045)
    h_El_Pi.SetLineColor(r.kBlue)
    h_El_Pi.Draw('histo')
    h_El_K = h_K.ProjectionY()
    h_El_K.SetLineColor(r.kRed)
    h_El_K.Draw('histo sames')
    c.cd(3)
    h_Mmiss_Pi = h_Pi.ProjectionZ()
    h_Mmiss_Pi.SetTitle("")
    h_Mmiss_Pi.SetTitle("")
    h_Mmiss_Pi.GetYaxis().SetRangeUser(0,6.5E3)
    h_Mmiss_Pi.GetYaxis().SetLabelSize(0.045)
    h_Mmiss_Pi.GetXaxis().SetLabelSize(0.045)
    h_Mmiss_Pi.GetXaxis().SetTitleSize(0.045)
    h_Mmiss_Pi.SetLineColor(r.kBlue)
    h_Mmiss_Pi.Draw('histo')
    h_Mmiss_K = h_K.ProjectionZ()
    h_Mmiss_K.SetLineColor(r.kRed)
    h_Mmiss_K.Draw('histo sames')
    c.SaveAs('plots/MISIDTemplates_'+sample+'.C')

    print('Data pi MISID Yield: {:.2f}'.format(h_Pi.Integral()))
    print('Data K MISID Yield: {:.2f}'.format(h_K.Integral()))
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
