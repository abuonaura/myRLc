import ROOT as r
import numpy as np
import os.path
import os,sys,getopt,time

datadir = '$FILEDIR/MISID/SameSign/'
sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

nbins = 20
h_nTracks_K = r.TH1F('h_nTracks_K',';nTracks;',nbins,0,500)
h_nTracks_Pi = r.TH1F('h_nTracks_Pi',';nTracks;',nbins,0,500)
h_p_K = r.TH1F('h_p_K',';p(MeV/c);',nbins,0,1.4E5)
h_p_Pi = r.TH1F('h_p_Pi',';p(MeV/c);',nbins,0,1.4E5)
h_eta_K = r.TH1F('h_eta_K',';#eta;',nbins,0,8)
h_eta_Pi = r.TH1F('h_eta_Pi',';#eta;',nbins,0,8)



def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h


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

    polarities = ['MagUp','MagDown']
    for polarity in polarities:
        fname_Ksample = datadir+'K_sample_'+polarity+sample_suffix[sample]
        fname_Pisample = datadir+'Pi_sample_'+polarity+sample_suffix[sample]
        f_Ksample = r.TFile(fname_Ksample,"READ")
        t_Ksample = f_Ksample.Get('DecayTree')
        f_Pisample = r.TFile(fname_Pisample,"READ")
        t_Pisample = f_Pisample.Get('DecayTree')
        
        for nevt in range(t_Ksample.GetEntries()):
            t_Ksample.GetEntry(nevt)
            nTracks = t_Ksample.nTracks
            sig_sw = t_Ksample.sw_sig
            p = t_Ksample.mu_P
            pl = t_Ksample.mu_PZ
            eta = (1./2.)*np.log((p+pl)/(p-pl))
            h_nTracks_K.Fill(nTracks,sig_sw)
            h_p_K.Fill(p,sig_sw)
            h_eta_K.Fill(eta,sig_sw)

        for nevt in range(t_Pisample.GetEntries()):
            t_Pisample.GetEntry(nevt)
            nTracks = t_Pisample.nTracks
            p = t_Pisample.mu_P
            pl = t_Pisample.mu_PZ
            eta = (1./2.)*np.log((p+pl)/(p-pl))
            sig_sw = t_Pisample.sw_sig
            h_nTracks_Pi.Fill(nTracks,sig_sw)
            h_p_Pi.Fill(p,sig_sw)
            h_eta_Pi.Fill(eta,sig_sw)

    h_nTracks_K = ScaleHisto(h_nTracks_K,1)
    h_nTracks_Pi = ScaleHisto(h_nTracks_Pi,1)
    h_p_K = ScaleHisto(h_p_K,1)
    h_p_Pi = ScaleHisto(h_p_Pi,1)
    h_eta_K = ScaleHisto(h_eta_K,1)
    h_etas_Pi = ScaleHisto(h_eta_Pi,1)

    c = r.TCanvas('c','',1500,500)
    c.Divide(3,1)
    c.cd(1)
    h_nTracks_K.SetMarkerStyle(20)
    h_nTracks_K.SetMarkerSize(0.7)
    h_nTracks_K.Draw('P0')
    h_nTracks_Pi.SetMarkerStyle(20)
    h_nTracks_Pi.SetMarkerSize(0.7)
    h_nTracks_Pi.SetMarkerColor(r.kRed)
    h_nTracks_Pi.Draw('P0 sames')
    l = r.TLegend(0.1,0.75,0.4,0.9)
    l.AddEntry(h_nTracks_K,"K_sample","p")
    l.AddEntry(h_nTracks_Pi,"Pi_sample","p")
    l.Draw()
    c.cd(2)
    h_p_K.SetMarkerStyle(20)
    h_p_K.SetMarkerSize(0.7)
    h_p_K.Draw('P0')
    h_p_Pi.SetMarkerStyle(20)
    h_p_Pi.SetMarkerSize(0.7)
    h_p_Pi.SetMarkerColor(r.kRed)
    h_p_Pi.Draw('P0 sames')
    c.cd(3)
    h_eta_K.SetMarkerStyle(20)
    h_eta_K.SetMarkerSize(0.7)
    h_eta_K.Draw('P0')
    h_eta_Pi.SetMarkerStyle(20)
    h_eta_Pi.SetMarkerSize(0.7)
    h_eta_Pi.SetMarkerColor(r.kRed)
    h_eta_Pi.Draw('P0 sames')
