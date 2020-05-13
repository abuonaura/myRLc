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




datadir = '$FILEDIR/'
polarities=['MagUp','MagDown']
particles=['K','Pi']

sample_suffix = {'full':'_sw.root', 'iso':'_iso_sw.root','Kenriched':'_Kenr_sw.root'}
sample_suffixCF = {'full':'_sw_withCF.root', 'iso':'_iso_sw_withCF.root','Kenriched':'_Kenr_sw_withCF.root'}
suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

ISOBDTcut = 0.35
ISOBDT2cut = 0.2

h_data = r.TH3F('h_data',';E_{#mu} (MeV);q^{2} (MeV^{2});Mmiss^{2} (MeV^{2}',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)

h_muMISID_Pi = r.TH3F('h_muMISID_Pi',';E_{#mu} (MeV);q^{2} (MeV^{2});Mmiss^{2} (MeV^{2}',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)

h_muMISID_K = r.TH3F('h_muMISID_K',';E_{#mu} (MeV);q^{2} (MeV^{2});Mmiss^{2} (MeV^{2}',20,0,2500,20,-2E6,18E6,20,-2E6,18E6)

h_muMISID = {'K':h_muMISID_K,'Pi':h_muMISID_Pi}

fitvar = ['FitVar_El_mLc','FitVar_q2_mLc','FitVar_Mmiss2_mLc']

def RoundHisto(h):
    for i in range(h.GetNbinsX()):
        for j in range(h.GetNbinsY()):
            for k in range(h.GetNbinsZ()):
                n = h.GetBinContent(i,j,k)
                if n<1:
                    nbin = h.GetBin(i,j,k)
                    h.SetBinContent(nbin, 0.)
    return h



def GetTemplatesData(sample):
    for polarity in polarities:
        if sample!='Kenriched':
            datafname = datadir+'Data/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        else:
            datafname = datadir+'ControlSamples/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        dataf = r.TFile(datafname, 'READ')
        datat = dataf.Get('DecayTree')
        for i in range(datat.GetEntries()):
            datat.GetEntry(i)
            h_data.Fill(datat.FitVar_El_mLc,datat.FitVar_q2_mLc,datat.FitVar_Mmiss2_mLc,datat.sw_sig)  
    return h_data

def GetTemplatesMISID(particle, sample):
    for polarity in polarities:
        misidfname = datadir+'/MISID/'+particle+'_sample_'+polarity+sample_suffixCF[sample]
        misidf = r.TFile(misidfname,'READ')
        misidt = misidf.Get('DecayTree')
        for i in range(misidt.GetEntries()):
            misidt.GetEntry(i)
            h_muMISID[particle].Fill(misidt.FitVar_El_mLc,misidt.FitVar_q2_mLc,misidt.FitVar_Mmiss2_mLc,misidt.sw_sig * misidt.w_recomu_CF *10)
    return h_muMISID[particle]

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

def ProduceTemplates(sample):
    h_data = GetTemplatesData(sample)
    for particle in particles:
        h_muMISID[particle] = GetTemplatesMISID(particle,sample)
    return h_data, h_muMISID

def SaveTemplates2file(filename,h_data, h_muMISID,sample):
    f = r.TFile(filename,'recreate')
    h_data, h_muMISID = ProduceTemplates(sample)
    h_data = RoundHisto(h_data)
    h_data.SetDirectory(f)
    for particle in particles:
        h_muMISID[particle] = RoundHisto(h_muMISID[particle])
        h_muMISID[particle].SetDirectory(f)
    f.Write()
    f.Close()

def ComputeMISIDFractions(filenameIN, filenameOUT):
    f = r.TFile(filenameIN,'read')
    f1 = r.TFile(filenameOUT,'recreate')
    h_data = f.Get('h_data')
    h_muMISID_Pi = f.Get('h_muMISID_Pi')
    h_muMISID_K = f.Get('h_muMISID_K')
    #h_fPi = fraction of pions misidentified as muons in each bin
    h_fPi = r.TH3F('h_fPi',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,
18E6,20,-2E6,18E6)
    h_fPi.Divide(h_muMISID_Pi, h_data)
    #h_fK = fraction of kaons misidentified as muons in each bin
    h_fK = r.TH3F('h_fK',';E_{#mu} (MeV); q^{2} (MeV^{2}); Mmiss^{2} (MeV^{2});',20,0,2500,20,-2E6,
18E6,20,-2E6,18E6)
    h_fK.Divide(h_muMISID_K, h_data)
    #save to file
    h_fPi.SetDirectory(f1)
    h_fK.SetDirectory(f1)
    f1.Write()
    f1.Close()
    f.Close()

def AddMISIDweights(ifile,ofile, fractionfile, sample):
    #---> Open histrograms with muon misid fractions
    f1 = r.TFile(fractionfile,'read')
    h_fPi = f1.Get('h_fPi')
    h_fK = f1.Get('h_fK')

    f = r.TFile(ifile,'READ')
    t = f.Get('DecayTree')

    of = r.TFile(ofile,'recreate')
    print ("Copying the original tree ...")
    ot = t.CloneTree()
    print ("Tree copied.")

    #Probability of muon misidentification
    w_MISID = np.zeros(1, dtype=float)
    b_MISID = ot.Branch('w_MISID',w_MISID,'w_MISID/D')

    #-----> Loop over the entries of the data sample
    for nevt in range(ot.GetEntries()):
        ot.GetEntry(nevt)
        E = ot.FitVar_El_mLc
        q2 = ot.FitVar_q2_mLc
        Mmiss = ot.FitVar_Mmiss2_mLc
        #-----> Look at the ID of the corresponding bin in the histos
        ibinK = h_fK.FindFixBin(E,q2,Mmiss)
        ibinPi = h_fPi.FindFixBin(E,q2,Mmiss)
        w_MISID[0] = (1 - h_fK.GetBinContent(ibinK) - h_fPi.GetBinContent(ibinPi))

        b_MISID.Fill()

    of.Write()
    of.Close()
    f.Close()
    f1.Close()


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


    filename = datadir+'MISID/IntermediateFiles_OS/Templates'+suffix[sample]
    fractionfile = datadir + 'MISID/IntermediateFiles_OS/muMISIDfractions'+suffix[sample]

    
    if os.path.isfile(filename) and restart==False:
        print('Template file already created')
    else:
        SaveTemplates2file(filename,h_data,h_muMISID,sample)
    print('Moving on to MISID fraction computation')
    ComputeMISIDFractions(filename, fractionfile)
    for polarity in polarities:
        print('>>>>>   Processing:          ', polarity )
        if sample!='Kenriched':
            outdatafile = datadir+'Data/Lb_Data_'+polarity+'_sw_noMISID'+suffix[sample]
            datafname = datadir+'Data/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        else:
            outdatafile = datadir+'ControlSamples/Lb_Data_'+polarity+'_sw_noMISID'+suffix[sample]
            datafname = datadir+'ControlSamples/Lb_Data_'+polarity+'_reduced_preselected'+sample_suffix[sample]
        print('- Reading data file: '+datafname)
        print('- Output data file: '+outdatafile)
        AddMISIDweights(datafname, outdatafile, fractionfile, sample)
            


