'''
Autor: Annarita Buonaura
Date: 30 December 2018

Description: Fills 3D histo with full Pi(K) samples and again reweighteng each evt for the ratio effMISID/effID and saves these histos to root file

'''

import ROOT as r
import numpy as np
import os.path
import os,sys,getopt,time

#INPUT OPTION DEFINITION
#full -> full sample, no cut on ISObdt
#iso -> isolated region: Lb_ISOLATION_BDT<0.35
#Kenriched -> anti-isolated region kaon enriched: (Lb_ISOLATION_BDT>0.35 && Lb_ISOLATION_PIDK>4.) || (Lb_ISOLATION_BDT2>0.35 && Lb_ISOLATION_PIDK2>4.)
#------
#6.06.2019
#Kenriched -> anti-isolated region kaon enriched: (Lb_ISOLATION_BDT>0.35 && Lb_ISOLATION_BDT2>0.35) &&(Lb_ISOLATION_PIDK>4.|| Lb_ISOLATION_PIDK2>4.)
#------
#4.07.2019 
#Modified anti-isolated region k enriched definition: (Lb_ISOLATION_BDT>"+str(ISOBDTcut)+"&& Lb_ISOLATION_BDT2>" +str(ISOBDT2cut)+")&&((Lb_ISOLATION_PIDK>4.&&(Lb_ISOLATION_CHARGE==mu_ID/13 ||(Lb_ISOLATION_CHARGE==-mu_ID/13 && Lb_ISOLATION_PIDp - Lb_ISOLATION_PIDK<0.))) || (Lb_ISOLATION_PIDK2>4.&&(Lb_ISOLATION_CHARGE2==mu_ID/13 ||(Lb_ISOLATION_CHARGE2==-mu_ID/13 && Lb_ISOLATION_PIDp2 - Lb_ISOLATION_PIDK2<0.))))


sample_suffix = {'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
suffix = {'iso':'_iso.root','Kenriched':'_Kenr.root'}

datadir = '/disk/lhcb_data2/RLcMuonic2016/'

polarities=['MagUp','MagDown']
particles=['K','Pi']


def CreateHisto(polarity,particle,sample):
    fname= datadir+'/MISID/SameSign/'+particle+'_sample_'+polarity+sample_suffix[sample]
    f = r.TFile(fname,'READ')
    t = f.Get('DecayTree')

    h = r.TH3F('h_'+particle+'_'+polarity,';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_misid = r.TH3F('h_misid'+particle+'_'+polarity,';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        h.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig)
        h_misid.Fill(t.FitVar_El_mLc,t.FitVar_q2_mLc, t.FitVar_Mmiss2_mLc, t.sw_sig*t.w_recowhad)
    
    h.SetDirectory(0)
    h_misid.SetDirectory(0)

    return h,h_misid

def SaveHistos2file(sample):
    print('Saving Histos 2 file!')
    f = r.TFile(datadir+'/MISID/IntermediateFiles_SS/HistosK2Pi'+suffix[sample],'recreate')
    misidfname = {'MagUp':{'K':'','Pi':''},'MagDown':{'K':'','Pi':''}}
    for particle in particles:
        for polarity in polarities:
            misidfname[polarity][particle] = datadir+'/MISID/SameSign/'+particle+'_sample_'+polarity+sample_suffix[sample]

    h_Pi_MagUp,h_misidPi_MagUp = CreateHisto('MagUp','Pi',sample)
    h_Pi_MagDown,h_misidPi_MagDown = CreateHisto('MagDown','Pi',sample)
    h_K_MagUp, h_misidK_MagUp = CreateHisto('MagUp','K',sample)
    h_K_MagDown, h_misidK_MagDown = CreateHisto('MagDown','K',sample)
    
    h_Pi = r.TH3F('h_Pi',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_K = r.TH3F('h_K',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_Pi.Add(h_Pi_MagUp, h_Pi_MagDown)
    h_Pi.SetDirectory(f)
    h_K.Add(h_K_MagUp, h_K_MagDown)
    h_K.SetDirectory(f)
    h_misidPi = r.TH3F('h_misidPi',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_misidK = r.TH3F('h_misidK',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_misidPi.Add(h_misidPi_MagUp, h_misidPi_MagDown)
    h_misidK.Add(h_misidK_MagUp, h_misidK_MagDown)
    h_misidPi.SetDirectory(f)
    h_misidK.SetDirectory(f)
    f.Write()
    f.Close()
    return

def ComputeKPiMisidFractions(sample):
    f = r.TFile(datadir+'/MISID/IntermediateFiles_SS/HistosK2Pi'+suffix[sample],'read')
    f1 = r.TFile(datadir+'/MISID/IntermediateFiles_SS/FractionsK2Pi'+suffix[sample],'recreate')
    
    h_Pi = f.Get('h_Pi')
    h_misidPi = f.Get('h_misidPi')
    h_K = f.Get('h_K')
    h_misidK = f.Get('h_misidK')
    #---> h_fPi = fraction of Pions misidentified as kaons in each bin
    h_fPi = r.TH3F('h_fPi',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_fPi.Divide(h_misidPi, h_K)
    #---> h_fK = fraction of Kaons misidentified as Pions in each bin
    h_fK = r.TH3F('h_fK',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)
    h_fK.Divide(h_misidK, h_Pi)
    h_fPi.SetDirectory(f1)
    h_fK.SetDirectory(f1)

    f1.Write()
    f1.Close()
    f.Close()

    
def AddMISIDweightsWCF(ifile,ofile, polarity, particle,sample):
    #---> Open histrograms with crossfeed fractions
    f1 = r.TFile(datadir+'/MISID/IntermediateFiles_SS/FractionsK2Pi'+suffix[sample],'read')

    h_fPi = f1.Get('h_fPi')
    h_fK = f1.Get('h_fK')

    f = r.TFile(ifile,'READ')
    t = f.Get('DecayTree')

    of = r.TFile(ofile,'recreate')
    print ("Copying the original tree ...")
    ot = t.CopyTree('')
    print ("Tree copied.")

    #Probability of reconstructing the particle as a muon when accounting for CrossFeed
    w_recomu_CF = np.zeros(1, dtype=float)
    b_recomu_CF = ot.Branch('w_recomu_CF',w_recomu_CF,'w_recomu_CF/D')

    #-----> Loop over the entries of the subsample
    for nevt in range(ot.GetEntries()):
        ot.GetEntry(nevt)
        E = ot.FitVar_El_mLc
        q2 = ot.FitVar_q2_mLc
        Mmiss = ot.FitVar_Mmiss2_mLc

        #-----> Look at the ID of the corresponding bin in the histos
        if(particle=='Pi'):
            ibin = h_fK.FindFixBin(E,q2,Mmiss)
            w_recomu_CF[0] = ot.w_recomu*(1-h_fK.GetBinContent(ibin))
        if(particle=='K'):
            ibin = h_fPi.FindFixBin(E,q2,Mmiss)
            w_recomu_CF[0] = ot.w_recomu*(1-h_fPi.GetBinContent(ibin))

        b_recomu_CF.Fill()
    
    of.Write()
    of.Close()
    f.Close()
    f1.Close()

if __name__ == '__main__':
    restart=False
    opts, args = getopt.getopt(sys.argv[1:], "",["iso","Kenriched","restart"])
    print (opts,args)
    for o, a in opts:
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'
        if o in ("--restart",):
            restart = True

    print('>>>   Evaluating cross-feed for sample: ', sample)
    
    if (os.path.isfile(datadir+'/MISID/IntermediateFiles_SS/HistosK2Pi'+suffix[sample]) and restart==False) or restart==False:
        print('File with histograms already exists!')
    else:
        SaveHistos2file(sample)
    ComputeKPiMisidFractions(sample)
    for polarity in polarities:
        for particle in particles:
            print('>>>>>   Processing:          ', particle, ' ', polarity )
            misidfname = datadir+'/MISID/SameSign/'+particle+'_sample_'+polarity+sample_suffix[sample]
            outfname = misidfname[0:-5]+'_withCF.root'
            AddMISIDweightsWCF(misidfname,outfname,polarity, particle, sample)
