'''
Author: Annarita Buonaura
Data: 5 April 2019
Purpose: Draw fit templates of dataSS and misidSS samples for the different isolation cases (no iso, iso applied, anti-isolated sample)
+ Draw fit templates of data and (dataSS) combinatorial bkg
'''

import ROOT as r
import os, sys, getopt

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'

#polarities=['MagUp','MagDown']
polarities=['MagDown']
particles=['K','Pi']
sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
sample_suffix_misid = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}
sample_suffix_comb = {'full':'_sw_noMISID.root', 'iso':'_sw_noMISID_iso.root','Kenriched':'_sw_noMISID_Kenr.root'}


def makeTemplateCombinatorial(polarity,sample):
    histo = r.TH3F('h_'+polarity+'Comb','; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    if sample!='Kenriched':
        datafname = datadir + 'Data/Lb_DataSS_'+polarity+sample_suffix_comb[sample]
    else:
        datafname = datadir + 'ControlSamples/Lb_DataSS_'+polarity+sample_suffix_comb[sample]
    dataf = r.TFile(datafname, 'Read')
    datat = dataf.Get('DecayTree')
    print(datat.GetEntries())
    for i in range(datat.GetEntries()):
        datat.GetEntry(i)
        if datat.Lb_M<5620:
            histo.Fill(datat.FitVar_q2_mLc/1E6,datat.FitVar_El_mLc, datat.FitVar_Mmiss2_mLc/1E6,datat.sw_sig*datat.w_MISID)
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
    
    #Plotting data and combinatorial bkg templates
    h_MagUpComb = makeTemplateCombinatorial('MagUp',sample)
    h_MagDownComb = makeTemplateCombinatorial('MagDown',sample)
    h_comb = r.TH3F('h_comb','3D histo Comb bkg; q^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_comb.Add(h_MagUpComb, h_MagDownComb)


    r.gStyle.SetOptStat(0)
    c1 = r.TCanvas('c1','Combinatorial BKg',1600,500)
    c1.Divide(3,1)
    c1.cd(1)
    h_q2_comb = h_comb.ProjectionX()
    h_q2_comb.SetTitle("")
    h_q2_comb.SetLineColor(r.kBlack)
    h_q2_comb.GetYaxis().SetRangeUser(0,2.2E4)
    h_q2_comb.GetYaxis().SetLabelSize(0.045)
    h_q2_comb.GetXaxis().SetLabelSize(0.045)
    h_q2_comb.GetXaxis().SetTitleSize(0.045)
    h_q2_comb.Draw('histo')
    c1.cd(2)
    h_El_comb = h_comb.ProjectionY()
    h_El_comb.SetTitle("")
    h_El_comb.SetLineColor(r.kBlack)
    h_El_comb.GetYaxis().SetRangeUser(0,1.E4)
    h_El_comb.GetYaxis().SetLabelSize(0.045)
    h_El_comb.GetXaxis().SetLabelSize(0.045)
    h_El_comb.GetXaxis().SetTitleSize(0.045)
    h_El_comb.Draw('histo')
    c1.cd(3)
    h_Mmiss_comb = h_comb.ProjectionZ()
    h_Mmiss_comb.SetTitle("")
    h_Mmiss_comb.SetLineColor(r.kBlack)
    h_Mmiss_comb.GetYaxis().SetRangeUser(0,7E3)
    h_Mmiss_comb.GetYaxis().SetLabelSize(0.045)
    h_Mmiss_comb.GetXaxis().SetLabelSize(0.045)
    h_Mmiss_comb.GetXaxis().SetTitleSize(0.045)
    h_Mmiss_comb.Draw('histo')
    c1.SaveAs('plots/TemplatesONLYComb_'+sample+'.C')
    

