'''
Autor: Annarita Buonaura
Date: 13 December 2018

Description: Fills 3D histo with full Pi(K) samples and again reweighteng each evt for the ratio effMISID/effID and saves these histos to root file

'''

import ROOT as r
import numpy as np
import os,sys,getopt,time

#INPUT OPTION DEFINITION
#full -> full sample, no cut on ISObdt
#iso -> isolated region: Lb_ISOLATION_BDT<0.35
#anti-iso -> anti-isolated region kaon enriched: (Lb_ISOLATION_BDT>0.35 && Lb_ISOLATION_PIDK>4.) || (Lb_ISOLATION_BDT2>0.35 && Lb_ISOLATION_PIDK2>4.)
#------
#6.06.2019
#anti-iso -> anti-isolated region kaon enriched: (Lb_ISOLATION_BDT>0.35 && Lb_ISOLATION_BDT2>0.35) &&(Lb_ISOLATION_PIDK>4.|| Lb_ISOLATION_PIDK2>4.)
#------
#4.07.2019 
#Modified anti-isolated region k enriched definition: (Lb_ISOLATION_BDT>"+str(ISOBDTcut)+"&& Lb_ISOLATION_BDT2>" +str(ISOBDT2cut)+")&&((Lb_ISOLATION_PIDK>4.&&(Lb_ISOLATION_CHARGE==mu_ID/13 ||(Lb_ISOLATION_CHARGE==-mu_ID/13 && Lb_ISOLATION_PIDp - Lb_ISOLATION_PIDK<0.))) || (Lb_ISOLATION_PIDK2>4.&&(Lb_ISOLATION_CHARGE2==mu_ID/13 ||(Lb_ISOLATION_CHARGE2==-mu_ID/13 && Lb_ISOLATION_PIDp2 - Lb_ISOLATION_PIDK2<0.))))


datadir = '$FILEDIR/'
polarities=['MagUp','MagDown']
particles=['K','Pi']

sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
sample_suffixCF = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}
suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

def ComputeMISIDfraction_NoCF(particle, sample):
    nmisid = {'MagUp':0.0,'MagDown':0.0}
    nmisidTOT =0.
    nall = {'MagUp':0.0, 'MagDown':0.0}
    nallTOT=0.

    for polarity in polarities:
        if sample!='Kenriched':
            datafname = datadir+'Data/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        else:
            datafname = datadir+'ControlSamples/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]

        dataf = r.TFile(datafname, 'READ')
        datat = dataf.Get('DecayTree')

        misidfname = datadir+'/MISID/OppositeSign/'+particle+'_sample_'+polarity+sample_suffix[sample]
        misidf = r.TFile(misidfname,'READ')
        misidt = misidf.Get('DecayTree')

        for i in range(misidt.GetEntries()):
            misidt.GetEntry(i)
            nmisid[polarity] = nmisid[polarity]+(misidt.sw_sig * misidt.w_recomu *10)

        for j in range(datat.GetEntries()):
            datat.GetEntry(j)
            nall[polarity] = nall[polarity]+datat.sw_sig


        nmisidTOT = nmisidTOT+nmisid[polarity]
        nallTOT = nallTOT+nall[polarity]

    print('')
    print('Number of sweighted data evts: {}'.format(nallTOT))
    print('Number of {} reco as muons: {}'.format(particle, nmisidTOT))
    print('Fraction of misidentified {} events : {:.4f}%'.format(particle, (nmisidTOT/nallTOT)*100))
    print('')



def ComputeCFfraction(particle, sample):
    f = r.TFile(datadir+'MISID/IntermediateFiles_OS/HistosK2Pi'+suffix[sample],'read')
    
    if(particle=='Pi'):
        h = f.Get('h_Pi')
        h_misid = f.Get('h_misidK')
    if(particle=='K'):
        h = f.Get('h_K')
        h_misid = f.Get('h_misidPi')

    nmisid=0.
    nall=0.

    for i in range(1, h_misid.GetNbinsX() + 1):
        for j in range(1, h_misid.GetNbinsY()+1):
            for k in range(1,h_misid.GetNbinsZ()+1):
                ibin = h_misid.GetBin(i,j,k)
                nmisid = nmisid+h_misid.GetBinContent(ibin)
    for i in range(1, h.GetNbinsX() + 1):
        for j in range(1, h.GetNbinsY()+1):
            for k in range(1,h.GetNbinsZ()+1):
                ibin = h.GetBin(i,j,k)
                nall = nall+h.GetBinContent(ibin)

    if particle=='Pi':
        print ('Number of kaons reconstructed as pions: {:.2f}'.format(nmisid))
        print ('Number of pions: {:.2f}'.format(nall))
    if particle=='K':
        print ('Number of pions reconstructed as kaons: {:.2f}'.format(nmisid))
        print ('Number of kaons: {:.2f}'.format(nall))

    print ('Fraction of misidentified events: {:.4f}'.format(nmisid*1./nall))




def ComputeMISIDfraction_WithCF(particle, sample):
    nmisid = {'MagUp':0.0,'MagDown':0.0}
    nmisidTOT =0.
    nall = {'MagUp':0.0, 'MagDown':0.0}
    nallTOT=0.

    for polarity in polarities:
        if sample!='Kenriched':
            datafname = datadir+'Data/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        else:
            datafname = datadir+'ControlSamples/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        dataf = r.TFile(datafname, 'READ')
        datat = dataf.Get('DecayTree')

        misidfname = datadir+'/MISID/OppositeSign/'+particle+'_sample_'+polarity+sample_suffixCF[sample]
        misidf = r.TFile(misidfname,'READ')

        misidt = misidf.Get('DecayTree')

        for i in range(misidt.GetEntries()):
            misidt.GetEntry(i)
            nmisid[polarity] = nmisid[polarity]+(misidt.sw_sig * misidt.w_recomu_CF *10)
        for j in range(datat.GetEntries()):
            datat.GetEntry(j)
            nall[polarity] = nall[polarity]+datat.sw_sig

        #print(nall[polarity])
        nmisidTOT = nmisidTOT+nmisid[polarity]
        nallTOT = nallTOT+nall[polarity]
    
    print('')
    print('Number of sweighted data evts: {}'.format(nallTOT))
    print('Number of {} reco as muons: {}'.format(particle, nmisidTOT))
    print('Fraction of misidentified {} events (Crossfeed considered): {:.4f}%'.format(particle, (nmisidTOT/nallTOT)*100))
    print('')


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

    print('>>>>   Processing '+sample+' sample')
    for particle in particles:
        ComputeMISIDfraction_NoCF(particle, sample)
        print('')
        ComputeMISIDfraction_WithCF(particle, sample)
        print('')
        ComputeCFfraction(particle, sample)

