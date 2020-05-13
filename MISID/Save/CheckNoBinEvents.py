'''
    Author: Annarita Buonaura
    Date: March 25,2019

    Description: Checks how many events have no corresponding kinematic bin in PIDCalib tables
'''
import ROOT as r
import numpy as np
from ROOT import TFile, TTree, TBranch, TTreeReader, TH1F, TCanvas, TMath, gSystem

applyiso=True

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'

pathID = datadir+'HistoPID/IDeff/'

pathMisid = {'Mu':datadir+'HistoPID/MISIDeff/','K': datadir+'HistoPID/K2Pi/','Pi':datadir+'HistoPID/Pi2K/'}

particles = ['K','Pi']

polarities = ['MagUp','MagDown']

files = {'MagUp': datadir+'Data/Lb_FakeMu_MagUp_reduced_preselected_sw.root',
    'MagDown': datadir+'Data/Lb_FakeMu_MagDown_reduced_preselected_sw.root'}

cuts = {'K':'mu_PIDK >4.0',
        'Pi':'mu_PIDK<2',
        'P':'mu_PIDp < 0.0 && (mu_MC15TuneV1_ProbNNp - mu_MC15TuneV1_ProbNNk)>0.',
        'Mu':'mu_PIDmu>-200'}
K2pi_MisId = 'DLLK<2.0'
Pi2K_MisId = 'DLLK>4.0'

def GetPIDhisto(PIDfilename,histoname):
    PIDfile = r.TFile(PIDfilename,'READ')
    h = PIDfile.Get(histoname)
    h.SetDirectory(0)
    return h

def CheckHistoRanges(polarity, particle):
    reco_had_had = {'K':'DLLK >4.0 &&IsMuon==0','Pi':'DLLK<2 && IsMuon==0'}
    reco_mu_had = {'K':'DLLK >4.0 &&IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0', 'Pi':'DLLK<2 && IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0'}
    reco_wrong_had = {'K':'IsMuon==0 && ' + K2pi_MisId, 'Pi': 'IsMuon==0 && '+Pi2K_MisId}
    #-----> Open PID histograms
    fhistname = 'PerfHists_'+particle+'_Turbo16_'+polarity+'_P_ETA_nTracks_Brunel.root'
    hnameID = particle + '_'+reco_had_had[particle]+'_All'
    hnameMuMisID = particle + '_'+reco_mu_had[particle]+'_All'
    hnameHadMisID = particle + '_'+reco_wrong_had[particle]+'_All'

    hID = GetPIDhisto(pathID+fhistname,hnameID)
    print('>>>    Histo P('+particle+'|'+particle+')')
    print('p range: [', hID.GetXaxis().GetXmin(),',' , hID.GetXaxis().GetXmax(),']')
    print('Eta range: [', hID.GetYaxis().GetXmin(),',',  hID.GetYaxis().GetXmax(),']')
    print('Ntracks range: [', hID.GetZaxis().GetXmin(),',',  hID.GetZaxis().GetXmax(),']')
    hMuMISID = GetPIDhisto(pathMisid['Mu']+fhistname,hnameMuMisID)
    print('>>>    Histo P('+particle+'|mu)')
    print('p range: [', hMuMISID.GetXaxis().GetXmin(),',',  hMuMISID.GetXaxis().GetXmax(),']')
    print('Eta range: [', hMuMISID.GetYaxis().GetXmin(),',',  hMuMISID.GetYaxis().GetXmax(),']')
    print('Ntracks range: [', hMuMISID.GetZaxis().GetXmin(),',', hMuMISID.GetZaxis().GetXmax(),']')
    hHadMISID = GetPIDhisto(pathMisid[particle]+fhistname,hnameHadMisID)
    if particle=='Pi':
        print('>>>    Histo P('+particle+'|K)')
    else:
        print('>>>    Histo P('+particle+'|Pi)')
    print('p range: [', hHadMISID.GetXaxis().GetXmin(),',',  hHadMISID.GetXaxis().GetXmax(),']')
    print('Eta range: [', hHadMISID.GetYaxis().GetXmin(),',',  hHadMISID.GetYaxis().GetXmax(),']')
    print('Ntracks range: [', hHadMISID.GetZaxis().GetXmin(),',',  hHadMISID.GetZaxis().GetXmax(),']')






def CheckNumber0Events(ifile, applyiso, polarity, particle):
    reco_had_had = {'K':'DLLK >4.0 &&IsMuon==0','Pi':'DLLK<2 && IsMuon==0'}
    reco_mu_had = {'K':'DLLK >4.0 &&IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0', 'Pi':'DLLK<2 && IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0'}
    reco_wrong_had = {'K':'IsMuon==0 && ' + K2pi_MisId, 'Pi': 'IsMuon==0 && '+Pi2K_MisId}
    ISOBDTcut = 0.35
    isocut = "Lb_ISOLATION_BDT<"+str(ISOBDTcut)

    #-----> Open PID histograms
    fhistname = 'PerfHists_'+particle+'_Turbo16_'+polarity+'_P_ETA_nTracks_Brunel.root'
    hnameID = particle + '_'+reco_had_had[particle]+'_All'
    hnameMuMisID = particle + '_'+reco_mu_had[particle]+'_All'
    hnameHadMisID = particle + '_'+reco_wrong_had[particle]+'_All'

    hID = GetPIDhisto(pathID+fhistname,hnameID)
    hMuMISID = GetPIDhisto(pathMisid['Mu']+fhistname,hnameMuMisID)
    hHadMISID = GetPIDhisto(pathMisid[particle]+fhistname,hnameHadMisID)

    f = r.TFile(ifile,'READ')
    t = f.Get('DecayTree')

    outrange_max = outrange_min= emptybinID = emptybinMUID = emptybinHMISID=0

    of = r.TFile('prova.root','recreate')
    print ("Copying the original tree ...")
    if applyiso:
        ot = t.CopyTree(cuts[particle]+"&&"+isocut)
    else:
        ot = t.CopyTree(cuts[particle])

    print ("Tree copied.")

    #-----> Loop over the entries of the subsample
    for nevt in range(ot.GetEntries()):
        ot.GetEntry(nevt)
        nTracks = ot.nTracks
        p = ot.mu_P
        pl = ot.mu_PZ
        eta = (1./2.)*np.log((p+pl)/(p-pl))

        #-----> Look at the ID of the corresponding bin in the histos
        if p>hID.GetXaxis().GetXmax() or eta>hID.GetYaxis().GetXmax() or nTracks>hID.GetZaxis().GetXmax():
            outrange_max+=1
        if p<hID.GetXaxis().GetXmin() or eta<hID.GetYaxis().GetXmin() or nTracks<hID.GetZaxis().GetXmin():
            outrange_min+=1
        if (p<=hID.GetXaxis().GetXmax() and eta<=hID.GetYaxis().GetXmax() and nTracks<=hID.GetZaxis().GetXmax()) and (p>=hID.GetXaxis().GetXmin() or eta>=hID.GetYaxis().GetXmin() or nTracks>=hID.GetZaxis().GetXmin()):
            ibinID = hID.FindFixBin(p,eta,nTracks)
            ibinMuMISID = hMuMISID.FindFixBin(p,eta,nTracks)
            ibinHadMISID = hHadMISID.FindFixBin(p,eta,nTracks)

            #-----> Retrieve the efficiency for the corresponding bin
            prob_recoh = hID.GetBinContent(ibinID)
            if(prob_recoh==0): 
                emptybinID+=1
            prob_recomu = hMuMISID.GetBinContent(ibinMuMISID)
            if(prob_recomu==0): 
                emptybinMUID+=1
            prob_recowh = hHadMISID.GetBinContent(ibinHadMISID)
            if(prob_recowh==0): 
                emptybinHMISID+=1




    print('Total number of events: ', ot.GetEntries())
    print('Number of events out of range in P, eta OR nTracks: ', outrange_max, ' (above)     ', outrange_min, '  (below)')
    print('Number of events with empty P('+particle+'|'+particle+'): ', emptybinID)
    print('Number of events with empty P('+particle+'|mu): ', emptybinMUID)
    if particle=='Pi':
        print('Number of events with empty P('+particle+'|K): ', emptybinHMISID)
    if particle=='K':
        print('Number of events with empty P('+particle+'|Pi): ', emptybinHMISID)



for polarity in polarities:
    for particle in particles:
        print('Apply Isolation Bdtcut: ', applyiso)
        print('>>>>>   Processing:          ', particle, ' ', polarity )
        CheckHistoRanges(polarity,particle)
        #CheckNumber0Events(files[polarity],applyiso,polarity, particle)
