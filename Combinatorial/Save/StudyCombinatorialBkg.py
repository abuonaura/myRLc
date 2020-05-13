'''
Author: A. Buonaura
Date: 14th June 2019
 - STEP1: compare Data & DataSS template samples for m>m(Lb) after discarding MISID evts (=> reweight for w_MISID)
'''

import ROOT as r
from ROOT import TFile, TTree, TH1F, TH3D, TCanvas
import os, sys, getopt, time


datadir = '$FILEDIR/'
datatypes = ['Data','DataSS']
polarities = ['MagUp','MagDown']
particles = ['K','Pi']

sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
sample_suffixCF = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}
suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

LbMass = 5620.2

def ReturnDatafname(type, polarity, sample):
    if sample!='Kenriched':
        datafiles = {'Data':{'MagUp': datadir + 'Data/Lb_Data_MagUp_sw_noMISID'+suffix[sample],'MagDown': datadir+'Data/Lb_Data_MagDown_sw_noMISID'+suffix[sample]},'DataSS':{'MagUp': datadir + 'Data/Lb_DataSS_MagUp_sw_noMISID'+suffix[sample],'MagDown': datadir+'Data/Lb_DataSS_MagDown_sw_noMISID'+suffix[sample]}}
    else:
        datafiles = {'Data':{'MagUp': datadir + 'ControlSamples/Lb_Data_MagUp_sw_noMISID'+suffix[sample],'MagDown': datadir+'ControlSamples/Lb_Data_MagDown_sw_noMISID'+suffix[sample]},'DataSS':{'MagUp': datadir + 'ControlSamples/Lb_DataSS_MagUp_sw_noMISID'+suffix[sample],'MagDown': datadir+'ControlSamples/Lb_DataSS_MagDown_sw_noMISID'+suffix[sample]}}
    fname = datafiles[type][polarity]
    return fname

def ReturnMISIDfname(polarity, particle, sample):
    misidfilesCF = {'MagUp':{'K':datadir+'/MISID/OppositeSign/K_sample_MagUp'+sample_suffixCF[sample],'Pi':datadir+'/MISID/OppositeSign/Pi_sample_MagUp'+sample_suffixCF[sample]},'MagDown':{'K':datadir+'/MISID/OppositeSign/K_sample_MagDown'+sample_suffixCF[sample],'Pi':datadir+'/MISID/OppositeSign/Pi_sample_MagDown'+sample_suffixCF[sample]}}
    fname = misidfilesCF[polarity][particle]
    return fname

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h


def SaveHistos2file(sample):
    fnew = r.TFile('$FILEDIR/CombinatorialBkg/NonCorrectedHistos'+suffix[sample],'recreate')
    print('Creating file NonCorrectedHistos'+suffix[sample])
    #Fill the histograms for both Data and DataSS in the 3 Lb_M regions (above, below, full)
    h_data_a = r.TH3F('h_data_above',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,12000,20,-20E6,20E6,20,-20E6,20E6)
    h_data_b = r.TH3F('h_data_below',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,5000,20,-2E6,20E6,20,-2E6,20E6)
    h_data_f = r.TH3F('h_data_full',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,12000,20,-20E6,20E6,20,-20E6,20E6)
    h_LbM_data_a = TH1F('h_LbM_data_a','m(#Lambda_{b}) Data above range',20,5000,10000)
    #h_LbM_data_b = TH1F('h_LbM_data_b','m(#Lambda_{b}) Data below range',50,0,6000)
    h_LbM_data_f = TH1F('h_LbM_data_f','m(#Lambda_{b}) Data full range',50,0,10000)

    h_dataSS_a = r.TH3F('h_dataSS_above',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,12000,20,-20E6,20E6,20,-20E6,20E6)
    h_dataSS_b = r.TH3F('h_dataSS_below',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,5000,20,-2E6,20E6,20,-2E6,20E6)
    h_dataSS_f = r.TH3F('h_dataSS_full',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,12000,20,-20E6,20E6,20,-20E6,20E6)
    h_LbM_dataSS_a = TH1F('h_LbM_dataSS_a','m(#Lambda_{b}) DataSS above range',20,5000,10000)
    #h_LbM_dataSS_b = TH1F('h_LbM_dataSS_b','m(#Lambda_{b}) DataSS below range',50,0,6000)
    h_LbM_dataSS_f = TH1F('h_LbM_dataSS_f','m(#Lambda_{b}) DataSS full range',50,0,10000)


    for dtype in datatypes:
        for polarity in polarities:
            fname = ReturnDatafname(dtype,polarity,sample)
            f= TFile(fname,'READ')
            t = f.Get('DecayTree')

            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if dtype=='Data':
                    h_LbM_data_f.Fill(t.Lb_M, t.sw_sig)
                    h_data_f.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig)
                    if t.Lb_M>=LbMass and sample!='Kenriched':
                        h_data_a.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_MISID)
                        h_LbM_data_a.Fill(t.Lb_M, t.sw_sig*t.w_MISID)
                    if t.Lb_M<LbMass:
                        h_data_b.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig)
                        #h_LbM_data_b.Fill(t.Lb_M, t.sw_sig)
                if dtype=='DataSS':
                    h_LbM_dataSS_f.Fill(t.Lb_M, t.sw_sig)
                    h_dataSS_f.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig)
                    if t.Lb_M>=LbMass and sample!='Kenriched':
                        h_dataSS_a.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_MISID)
                        h_LbM_dataSS_a.Fill(t.Lb_M, t.sw_sig*t.w_MISID)
                    if t.Lb_M<LbMass:
                        h_dataSS_b.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_MISID)
                        #h_LbM_dataSS_b.Fill(t.Lb_M, t.sw_sig)
   
    if sample!='Kenriched':
        h_data_a.SetDirectory(fnew)
        h_LbM_data_a.SetDirectory(fnew)
    #h_LbM_data_b.SetDirectory(fnew)
    h_data_b.SetDirectory(fnew)
    h_data_f.SetDirectory(fnew)
    h_LbM_data_f.SetDirectory(fnew)

    if sample!='Kenriched':
        h_dataSS_a.SetDirectory(fnew)
        h_LbM_dataSS_a.SetDirectory(fnew)
    
    h_dataSS_b.SetDirectory(fnew)
    h_dataSS_f.SetDirectory(fnew)
    #h_LbM_dataSS_b.SetDirectory(fnew)
    h_LbM_dataSS_f.SetDirectory(fnew)

    fnew.Write()
    fnew.Close()
    return

if __name__ == '__main__':
    restart = False
    opts, args = getopt.getopt(sys.argv[1:], "",["full","iso","Kenriched","restart"])
    print (opts,args)
    for o, a in opts:
        if o in ("--full",):
            sample = 'full'
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'
        if o in ("--restart",):
            restart = True

    print(sample)

    print('')
    print('Processing '+sample+' sample')
    print('')

    #Save histos to a file (to avoid loosing too much time)
    if os.path.isfile('$FILEDIR/CombinatorialBkg/NonCorrectedHistos'+suffix[sample]) and restart==False:
        print('File with histograms already exists!')    
    else:
        SaveHistos2file(sample)

    
    fHistos = TFile('$FILEDIR/CombinatorialBkg/NonCorrectedHistos'+suffix[sample],'READ')
    if sample!='Kenriched':
        h_data_a = fHistos.Get('h_data_above')
    h_data_b = fHistos.Get('h_data_below')
    h_data_f = fHistos.Get('h_data_full')
    h_LbM_data_a = fHistos.Get('h_LbM_data_a')
    #h_LbM_data_b = fHistos.Get('h_LbM_data_b')
    h_LbM_data_f = fHistos.Get('h_LbM_data_f')

    if sample!='Kenriched':
        h_dataSS_a = fHistos.Get('h_dataSS_above')
    h_dataSS_b = fHistos.Get('h_dataSS_below')
    h_dataSS_f = fHistos.Get('h_dataSS_full')
    h_LbM_dataSS_a = fHistos.Get('h_LbM_dataSS_a')
    #h_LbM_dataSS_b = fHistos.Get('h_LbM_dataSS_b')
    h_LbM_dataSS_f = fHistos.Get('h_LbM_dataSS_f')

    #Project the Histos

    if sample!='Kenriched':
        h_El_data_a = h_data_a.ProjectionX()
        h_q2_data_a = h_data_a.ProjectionY()
        h_mmiss_data_a = h_data_a.ProjectionZ()
    
    h_El_data_b = h_data_b.ProjectionX()
    h_q2_data_b = h_data_b.ProjectionY()
    h_mmiss_data_b = h_data_b.ProjectionZ()

    h_El_data_f = h_data_f.ProjectionX()
    h_q2_data_f = h_data_f.ProjectionY()
    h_mmiss_data_f = h_data_f.ProjectionZ()
    
    if sample!='Kenriched':
        h_El_dataSS_a = h_dataSS_a.ProjectionX()
        h_q2_dataSS_a = h_dataSS_a.ProjectionY()
        h_mmiss_dataSS_a = h_dataSS_a.ProjectionZ()
    
    h_El_dataSS_b = h_dataSS_b.ProjectionX()
    h_q2_dataSS_b = h_dataSS_b.ProjectionY()
    h_mmiss_dataSS_b = h_dataSS_b.ProjectionZ()

    h_El_dataSS_f = h_dataSS_f.ProjectionX()
    h_q2_dataSS_f = h_dataSS_f.ProjectionY()
    h_mmiss_dataSS_f = h_dataSS_f.ProjectionZ()
    
    #Scale the histograms to 1
    if sample!='Kenriched':
        h_El_data_a = ScaleHisto(h_El_data_a,1)
        h_q2_data_a = ScaleHisto(h_q2_data_a,1)
        h_mmiss_data_a = ScaleHisto(h_mmiss_data_a,1)
        h_LbM_data_a = ScaleHisto(h_LbM_data_a,1)
    
    h_El_data_b = ScaleHisto(h_El_data_b,1)
    h_q2_data_b = ScaleHisto(h_q2_data_b,1)
    h_mmiss_data_b = ScaleHisto(h_mmiss_data_b,1)

    h_El_data_f = ScaleHisto(h_El_data_f,1)
    h_q2_data_f = ScaleHisto(h_q2_data_f,1)
    h_mmiss_data_f = ScaleHisto(h_mmiss_data_f,1)
    h_LbM_data_f = ScaleHisto(h_LbM_data_f,1)

    if sample!='Kenriched':
        h_El_dataSS_a = ScaleHisto(h_El_dataSS_a,1)
        h_q2_dataSS_a = ScaleHisto(h_q2_dataSS_a,1)
        h_mmiss_dataSS_a = ScaleHisto(h_mmiss_dataSS_a,1)
        h_LbM_dataSS_a = ScaleHisto(h_LbM_dataSS_a,1)
    
    h_El_dataSS_b = ScaleHisto(h_El_dataSS_b,1)
    h_q2_dataSS_b = ScaleHisto(h_q2_dataSS_b,1)
    h_mmiss_dataSS_b = ScaleHisto(h_mmiss_dataSS_b,1)

    h_El_dataSS_f = ScaleHisto(h_El_dataSS_f,1)
    h_q2_dataSS_f = ScaleHisto(h_q2_dataSS_f,1)
    h_mmiss_dataSS_f = ScaleHisto(h_mmiss_dataSS_f,1)
    h_LbM_dataSS_f = ScaleHisto(h_LbM_dataSS_f,1)

    #plot the full range kinematical templates superimposed
    c0 = TCanvas('c0','full range',1500,500)
    c0.Divide(3,1)
    c0.cd(1)
    h_El_data_f.SetMarkerStyle(20)
    h_El_data_f.SetMarkerSize(1)
    h_El_data_f.Draw('PE0')
    h_El_dataSS_f.SetLineColor(r.kRed)
    h_El_dataSS_f.Draw('hist sames')
    legend = r.TLegend(0.6,0.75,0.9,0.9)
    legend.AddEntry(h_El_dataSS_f,'Data SS', "l")
    legend.AddEntry(h_El_data_f,'Data OS', "pl")
    legend.Draw()
    c0.cd(2)
    h_q2_data_f.SetMarkerStyle(20)
    h_q2_data_f.SetMarkerSize(1)
    h_q2_data_f.Draw('PE0')
    h_q2_dataSS_f.SetLineColor(r.kRed)
    h_q2_dataSS_f.Draw('hist sames')
    c0.cd(3)
    h_mmiss_data_f.SetMarkerStyle(20)
    h_mmiss_data_f.SetMarkerSize(1)
    h_mmiss_data_f.Draw('PE0')
    h_mmiss_dataSS_f.SetLineColor(r.kRed)
    h_mmiss_dataSS_f.Draw('hist sames')

    #plot the full range kinematical templates superimposed
    if sample!='Kenriched':
        c0a = TCanvas('c0a','Compare above LbM range',1500,500)
        c0a.Divide(3,1)
        c0a.cd(1)
        h_El_data_a.SetMarkerStyle(20)
        h_El_data_a.SetMarkerSize(1)
        h_El_data_a.Draw('PE0')
        h_El_dataSS_a.SetLineColor(r.kRed)
        h_El_dataSS_a.Draw('hist sames')
        legend1 = r.TLegend(0.6,0.75,0.9,0.9)
        legend1.AddEntry(h_El_dataSS_a,'Data SS', "l")
        legend1.AddEntry(h_El_data_a,'Data OS', "pl")
        legend1.Draw()
        c0a.cd(2)
        h_q2_data_a.SetMarkerStyle(20)
        h_q2_data_a.SetMarkerSize(1)
        h_q2_data_a.Draw('PE0')
        h_q2_dataSS_a.SetLineColor(r.kRed)
        h_q2_dataSS_a.Draw('hist sames')
        c0a.cd(3)
        h_mmiss_data_a.SetMarkerStyle(20)
        h_mmiss_data_a.SetMarkerSize(1)
        h_mmiss_data_a.Draw('PE0')
        h_mmiss_dataSS_a.SetLineColor(r.kRed)
        h_mmiss_dataSS_a.Draw('hist sames')

    c1 = TCanvas('c1','Lb_M full range',1000,500)
    c1.Divide(2,1)
    c1.cd(1)
    h_LbM_data_f.Draw('hist')
    c1.cd(2)
    h_LbM_dataSS_f.SetLineColor(r.kRed)
    h_LbM_dataSS_f.Draw('hist')

    #Perform the ratio between the yields of histos above Lb_M
    h_LbMratio = r.TH1F('h_LbMratio',';m(#Lambda_{b}) (GeV^{2}/c^{4}); OS/SS',20,5000,10000)
    h_LbMratio.Divide(h_LbM_data_a,h_LbM_dataSS_a)
    if sample=='Kenriched':
        fit = r.TF1('fit','[0]*x+[1]',5000,10000)
    else:
        fit = r.TF1('fit','[0]*exp([1]*x+[2])',0,10000)

    if sample!='Kenriched':
        c2 = TCanvas('c2','Correction above Lb_M',1000,500)
        c2.Divide(2,1)
        c2.cd(1)
        h_LbM_dataSS_a.SetLineColor(r.kRed)
        h_LbM_dataSS_a.Draw('hist')
        h_LbM_data_a.Draw('hist sames')
        c2.cd(2)
        h_LbMratio.Draw()
        h_LbMratio.Fit(fit,"R+")
        fit.Draw('sames')
        a = fit.GetParameter(0)
        err_a = fit.GetParameter(0)
        b = fit.GetParameter(1)
        err_b = fit.GetParameter(1)
        c = fit.GetParameter(2)
        err_c = fit.GetParameter(2)
        #fit.Draw('same')
        print('Ratio: {:.4f} * exp({:.4f}*Lb_M + {:.4f})'.format(a,b,c))

        #Apply corrections to DataSS
        h_dataSS_a_c = r.TH3F('h_data_above_c',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,12000,20,-20E6,20E6,20,-20E6,20E6)
        h_dataSS_b_c = r.TH3F('h_data_below_c',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,5000,20,-2E6,20E6,20,-2E6,20E6)
                
    CombYield = 0
    CombYield_corr=0

    for polarity in polarities:
        fname = ReturnDatafname('DataSS', polarity,sample)
        f= TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_M>=LbMass:
                if sample!='Kenriched':
                    w_corr = fit.Eval(t.Lb_M)
                    h_dataSS_a_c.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_MISID*w_corr)
            if t.Lb_M<LbMass:
                if sample!='Kenriched':
                    w_corr = fit.Eval(t.Lb_M)
                    h_dataSS_b_c.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_MISID*w_corr)
                    CombYield_corr += t.sw_sig*t.w_MISID*w_corr
                CombYield += t.sw_sig*t.w_MISID

    #Evaluate signal yield:
    SigYield = 0
    for polarity in polarities:
        fname = ReturnDatafname('Data', polarity,sample)
        f= TFile(fname,'READ')
        t = f.Get('DecayTree')
        SigYield += float(t.Draw('Lb_M','sw_sig*(Lb_M<'+str(LbMass)+')','goff'))
    
    #Print Signal Yield
    print(' ******** YIELDS ******** ')
    print('Signal Yield: ', SigYield)
    print('Combinatorial Bkg Yield: ', CombYield)
    if sample!='Kenriched':
        print('Combinatorial Bkg (after r correction) Yield: ', CombYield_corr)

    #project the histos
    if sample!='Kenriched':
        h_El_dataSS_a_c = h_dataSS_a_c.ProjectionX()
        h_q2_dataSS_a_c = h_dataSS_a_c.ProjectionY()
        h_mmiss_dataSS_a_c = h_dataSS_a_c.ProjectionZ()
    
        h_El_dataSS_b_c = h_dataSS_b_c.ProjectionX()
        h_q2_dataSS_b_c = h_dataSS_b_c.ProjectionY()
        h_mmiss_dataSS_b_c = h_dataSS_b_c.ProjectionZ()

    #Scale them
        h_El_dataSS_a_c = ScaleHisto(h_El_dataSS_a_c,1)
        h_q2_dataSS_a_c = ScaleHisto(h_q2_dataSS_a_c,1)
        h_mmiss_dataSS_a_c = ScaleHisto(h_mmiss_dataSS_a_c,1)

        h_El_dataSS_b_c = ScaleHisto(h_El_dataSS_b_c,1)
        h_q2_dataSS_b_c = ScaleHisto(h_q2_dataSS_b_c,1)
        h_mmiss_dataSS_b_c = ScaleHisto(h_mmiss_dataSS_b_c,1)

    #plot Non corrected and corrected templates for DataSS comparing them to DataOS in both regions
    if sample!='Kenriched':
        c3 = TCanvas('c3','Checks correction above LbM',1500,500)
        c3.Divide(3,1)
        c3.cd(1)
        h_El_data_a.SetMarkerStyle(20)
        h_El_data_a.SetMarkerSize(1)
        h_El_data_a.Draw('PE0')
        h_El_dataSS_a.SetLineColor(r.kRed)
        h_El_dataSS_a.Draw('hist sames')
        h_El_dataSS_a_c.SetLineColor(r.kBlack)
        h_El_dataSS_a_c.Draw('hist sames')
        legend2 = r.TLegend(0.6,0.75,0.9,0.9)
        legend2.AddEntry(h_El_data_a,'Data OS', "pl")
        legend2.AddEntry(h_El_dataSS_a,'Data SS', "l")
        legend2.AddEntry(h_El_dataSS_a_c,'Data SS corr', "l")
        legend2.Draw()
        c3.cd(2)
        h_q2_data_a.SetMarkerStyle(20)
        h_q2_data_a.SetMarkerSize(1)
        h_q2_data_a.Draw('PE0')
        h_q2_dataSS_a.SetLineColor(r.kRed)
        h_q2_dataSS_a.Draw('hist sames')
        h_q2_dataSS_a_c.SetLineColor(r.kBlack)
        h_q2_dataSS_a_c.Draw('hist sames')
        c3.cd(3)
        h_mmiss_data_a.SetMarkerStyle(20)
        h_mmiss_data_a.SetMarkerSize(1)
        h_mmiss_data_a.Draw('PE0')
        h_mmiss_dataSS_a.SetLineColor(r.kRed)
        h_mmiss_dataSS_a.Draw('hist sames')
        h_mmiss_dataSS_a_c.SetLineColor(r.kBlack)
        h_mmiss_dataSS_a_c.Draw('hist sames')

        c4 = TCanvas('c4','Checks correction below LbM',1500,500)
        c4.Divide(3,1)
        c4.cd(1)
        h_El_data_b.SetMarkerStyle(20)
        h_El_data_b.SetMarkerSize(1)
        h_El_data_b.Draw('PE0')
        h_El_dataSS_b.SetLineColor(r.kRed)
        h_El_dataSS_b.Draw('hist sames')
        h_El_dataSS_b_c.SetLineColor(r.kBlack)
        h_El_dataSS_b_c.Draw('hist sames')
        legend3 = r.TLegend(0.6,0.75,0.9,0.9)
        legend3.AddEntry(h_El_data_a,'Data OS', "pl")
        legend3.AddEntry(h_El_dataSS_a,'Data SS', "l")
        legend3.AddEntry(h_El_dataSS_a_c,'Data SS corr', "l")
        legend3.Draw()
        c4.cd(2)
        h_q2_data_b.SetMarkerStyle(20)
        h_q2_data_b.SetMarkerSize(1)
        h_q2_data_b.Draw('PE0')
        h_q2_dataSS_b.SetLineColor(r.kRed)
        h_q2_dataSS_b.Draw('hist sames')
        h_q2_dataSS_b_c.SetLineColor(r.kBlack)
        h_q2_dataSS_b_c.Draw('hist sames')
        c4.cd(3)
        h_mmiss_data_b.SetMarkerStyle(20)
        h_mmiss_data_b.SetMarkerSize(1)
        h_mmiss_data_b.Draw('PE0')
        h_mmiss_dataSS_b.SetLineColor(r.kRed)
        h_mmiss_dataSS_b.Draw('hist sames')
        h_mmiss_dataSS_b_c.SetLineColor(r.kBlack)
        h_mmiss_dataSS_b_c.Draw('hist sames')

    else:
        c4 = TCanvas('c4','Checks correction below LbM',1500,500)
        c4.Divide(3,1)
        c4.cd(1)
        h_El_data_b.SetMarkerStyle(20)
        h_El_data_b.SetMarkerSize(1)
        h_El_data_b.Draw('PE0')
        h_El_dataSS_b.SetLineColor(r.kRed)
        h_El_dataSS_b.Draw('hist sames')
        legend3 = r.TLegend(0.6,0.75,0.9,0.9)
        legend3.AddEntry(h_El_data_b,'Data OS', "pl")
        legend3.AddEntry(h_El_dataSS_b,'Data SS', "l")
        legend3.Draw()
        c4.cd(2)
        h_q2_data_b.SetMarkerStyle(20)
        h_q2_data_b.SetMarkerSize(1)
        h_q2_data_b.Draw('PE0')
        h_q2_dataSS_b.SetLineColor(r.kRed)
        h_q2_dataSS_b.Draw('hist sames')
        c4.cd(3)
        h_mmiss_data_b.SetMarkerStyle(20)
        h_mmiss_data_b.SetMarkerSize(1)
        h_mmiss_data_b.Draw('PE0')
        h_mmiss_dataSS_b.SetLineColor(r.kRed)
        h_mmiss_dataSS_b.Draw('hist sames')
