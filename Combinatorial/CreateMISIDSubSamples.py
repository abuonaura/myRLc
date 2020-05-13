'''
    Author: Annarita Buonaura
    Date: May 8th, 2020

    Description: Creates MISID templates for same sign samples as needed for combinatorial bkg studies
'''
import ROOT as r
import numpy as np
from ROOT import TFile, TTree, TBranch, TTreeReader, TH1F, TCanvas, TMath, gSystem
gSystem.Load('libRooFit')
from ROOT import RooStats
from ROOT import RooFit as RF
from ROOT import RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData
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

suffix = {'iso':'_iso.root','Kenriched':'_Kenr.root'}

datadir = '/disk/lhcb_data2/RLcMuonic2016/'

pathID = datadir+'HistoPID/IDeff_Ghost/'
pathMisid = {'Mu':datadir+'HistoPID/MISIDeff_Ghost/','K': datadir+'HistoPID/K2Pi/','Pi':datadir+'HistoPID/Pi2K/'}

particles = ['K','Pi']

PIDcuts = {'K':'mu_PIDK >4.0 && mu_ProbNNghost < 0.2',
        'Pi':'mu_PIDK<2 && mu_ProbNNghost < 0.2'}

K2pi_MisId = 'DLLK<2.0'
Pi2K_MisId = 'DLLK>4.0'

ISOBDTcut = 0.35
ISOBDT2cut = 0.2



def GetDataset(t, *var):
    ds = r.RooDataSet("ds","ds",t,r.RooArgSet(*var))
    return ds

def addSignalBkgYields(start_sig,start_bkg,range_sig=[],range_bkg=[]):
    yields= r.RooArgList()
    sig = r.RooRealVar("nsig","N signal evts",start_sig,range_sig[0],range_sig[1])
    bkg = r.RooRealVar("nbkg","N bkg evts",start_bkg,range_bkg[0],range_bkg[1]);
    yields.add(sig)
    yields.add(bkg)
    return yields

def CreateMassFitModel(var, nsig,nbkg):
    w = r.RooWorkspace()

    #Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
    mean = r.RooRealVar("mean","mean of gaussians",2290,2270,2310) ;
    sigma1 = r.RooRealVar("sigma1","width of gaussians",3,0,10) ;
    sigma2 = r.RooRealVar("sigma2","width of gaussians",15,0,30) ;
    sig1 = r.RooGaussian("sig1","Signal component 1",var,mean,sigma1)
    sig2 = r.RooGaussian("sig2","Signal component 2",var,mean,sigma2)

    #Build exponential PDF
    l = r.RooRealVar("l", "slope of expo",-0.1, -1., 0.);
    bkg = r.RooExponential("bkg", " bkg with exponential PDF",var,l)

    #Sum the signal components into a composite signal p.d.f.
    sig1frac = r.RooRealVar("sig1frac","fraction of component 1 in signal",0.8,0.,1.)
    sig = r.RooAddPdf("sig","Signal",sig1,sig2,sig1frac)

    yields = r.RooArgList(nsig,nbkg)
    functions = r.RooArgList(sig,bkg)
    #Sum the composite signal and background
    model= r.RooAddPdf("model","g1+g2+a",functions,yields)

    getattr(w,'import')(model)
    return w

def GetWeightedDS(ds):
    ds_w = r.RooDataSet()
    ds_w = ds
    ds_w.SetName('dataWithSWeights')
    return ds_w

def PlotMassFit(model, ds, Lc_M,particle,polarity,sample):
    c = r.TCanvas()
    frame = Lc_M.frame(RF.Title('#Lambda_{c} mass peak distribution'))
    ds.plotOn(frame)
    model.plotOn(frame)
    model.paramOn(frame, RF.Layout(0.6,0.9,0.85))
    frame.Draw()
    c.Draw()
    imgsuff = suffix[sample][0:-5]+'.png'
    c.SaveAs('plots/Fit_'+particle+'_'+polarity+imgsuff)
    return

def AddSweights(ifile,polarity,particle,sample):
    #----> Open subsample files
    f = r.TFile(ifile,'read')
    t = f.Get('DecayTree')
    t.SetBranchStatus('sw_sig',0)
    t.SetBranchStatus('sw_bkg',0)

    #----> Create subsample files which will contain sweights
    of = r.TFile(ifile[0:-5]+'_sw.root','recreate')
    ot = t.CloneTree(0)
    ot.SetName('DecayTree')

    LcM_range = [2230,2330]
    Lc_M = r.RooRealVar('Lc_M','#Lambda_{c} mass',LcM_range[0],LcM_range[1])
    ds = GetDataset(t,Lc_M)

    #Create a workspace named 'w' with the model to fit the Lc_Mass peak
    nsig = r.RooRealVar("nsig","N signal evts",100000,0,500000)
    nbkg = r.RooRealVar("nbkg","N bkg evts",400000,0,500000)

    w = CreateMassFitModel(Lc_M,nsig,nbkg)
    model = w.pdf("model")
    result = model.fitTo(ds,RF.Extended(True),RF.Save())

    nsig = w.var('nsig')
    nbkg = w.var('nbkg')

    sData =r.RooStats.SPlot("sData","An SPlot",ds, model, r.RooArgList(nsig, nbkg))

    ds_w = GetWeightedDS(ds)

    #----> Create branches for sweights
    sw_bkg = np.zeros(1, dtype=float)
    sw_sig = np.zeros(1, dtype=float)

    ot.Branch("sw_bkg",sw_bkg,"sw_bkg/D")
    ot.Branch("sw_sig",sw_sig,"sw_sig/D")

    for i in range(ds_w.numEntries()):
        t.GetEntry(i)
        sw_bkg[0]= ds_w.get(i).getRealValue("nbkg_sw")
        sw_sig[0] = ds_w.get(i).getRealValue("nsig_sw")
        ot.Fill()

    ot.Write()
    of.Close()
    f.Close()

 #Save fit results to file
    of_fit = r.TFile(datadir+'/MISID/IntermediateFiles_SS/FitResults/FitResults_'+particle+'_'+polarity+ suffix[sample],'Recreate')
    result.Write()
    of_fit.Close()

    PlotMassFit(model, ds, Lc_M,particle,polarity,sample)
    
    return

def GetPIDhisto(PIDfilename,histoname):
    PIDfile = r.TFile(PIDfilename,'READ')
    h = PIDfile.Get(histoname)
    h.SetDirectory(0)
    return h

def AddMISIDweightNoCF(ifile, ofile, applyiso, polarity, particle):
    reco_had_had = {'K':'DLLK >4.0 &&IsMuon==0 && MC15TuneV1_ProbNNghost < 0.2','Pi':'DLLK<2 && IsMuon==0 && MC15TuneV1_ProbNNghost < 0.2'}
    reco_mu_had = {'K':'IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0 && DLLe<1 && MC15TuneV1_ProbNNghost < 0.2', 'Pi':'IsMuon==1 && DLLmu>2.0 && DLLmu-DLLK>2.0 && DLLmu-DLLp>2.0 && DLLe<1 && MC15TuneV1_ProbNNghost < 0.2'}
    reco_wrong_had = {'K':'IsMuon==0 && ' + K2pi_MisId, 'Pi': 'IsMuon==0 && '+Pi2K_MisId}

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

    of = r.TFile(ofile,'recreate')
    cut = PIDcuts[particle]
    print ("Copying the original tree ...")
    print('Applying the following cuts: '+cut)
    ot = t.CopyTree(cut)
    print ("Tree copied.")

    #Probability of reconstructing the particle as a muon
    w_recomu = np.zeros(1, dtype=float)
    b_recomu = ot.Branch('w_recomu',w_recomu,'w_recomu/D')
    
    #Probability of reconstructing the particle as a different hadron
    w_recowhad = np.zeros(1, dtype=float)
    b_recowhad = ot.Branch('w_recowhad',w_recowhad,'w_recowhad/D')

    #-----> Loop over the entries of the subsample
    for nevt in range(ot.GetEntries()):
        ot.GetEntry(nevt)
        nTracks = ot.nTracks
        p = ot.mu_P
        pl = ot.mu_PZ
        eta = (1./2.)*np.log((p+pl)/(p-pl))

        if p>hID.GetXaxis().GetXmax() or eta>hID.GetYaxis().GetXmax() or nTracks>hID.GetZaxis().GetXmax():
            w_recomu[0] = 0.0
            w_recowhad[0] = 0.0
            print('WARNING Event out of Hist axes range (P,eta,nTracks): {:.2f}, {:.2f}, {:.2f}'.format(p,eta,nTracks))
        if p<hID.GetXaxis().GetXmin() or eta<hID.GetYaxis().GetXmin() or nTracks<hID.GetZaxis().GetXmin():
            w_recomu[0] = 0.0
            w_recowhad[0] = 0.0
            print('WARNING Event out of Hist axes range (P,eta,nTracks): {:.2f}, {:.2f}, {:.2f}'.format(p,eta,nTracks))
        if (p<=hID.GetXaxis().GetXmax() and eta<=hID.GetYaxis().GetXmax() and nTracks<=hID.GetZaxis().GetXmax()) and (p>=hID.GetXaxis().GetXmin() or eta>=hID.GetYaxis().GetXmin() or nTracks>=hID.GetZaxis().GetXmin()):
            ibinID = hID.FindFixBin(p,eta,nTracks)
            ibinMuMISID = hMuMISID.FindFixBin(p,eta,nTracks)
            ibinHadMISID = hHadMISID.FindFixBin(p,eta,nTracks)

            #-----> Retrieve the efficiency for the corresponding bin
            prob_recoh = hID.GetBinContent(ibinID)
            prob_recomu = hMuMISID.GetBinContent(ibinMuMISID)
            prob_recowh = hHadMISID.GetBinContent(ibinHadMISID)

            if(prob_recoh!=0):
                w_recomu[0] = (prob_recomu/prob_recoh)
                w_recowhad[0] = (prob_recowh/prob_recoh)
            else:
                print ('WARNING in Hist p('+particle+'|'+particle+') Empty PIDCalib bin (P, eta, nTracks): {:.2f}, {:.2f}, {:.2f}'.format(p,eta,nTracks))
                w_recomu[0] = 0.0
                w_recowhad[0] = 0.0

        #-----> Fill the branches
        b_recomu.Fill()
        b_recowhad.Fill()

    #-----> Save modified subsample files & close
    ot.Write()
    of.Close()

    #-----> Add the sweights
    AddSweights(ofile,polarity,particle,sample)           

if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], "",["iso","Kenriched"])
    print (opts,args)
    for o, a in opts:
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'

    print (sample)
    polarities = ['MagUp','MagDown']
    files = {'MagUp': datadir+'Data/Lb_FakeMuSS_MagUp_preselected'+suffix[sample][0:-5]+'_sw.root',
            'MagDown': datadir+'Data/Lb_FakeMuSS_MagDown_preselected'+suffix[sample][0:-5]+'_sw.root'}

    for polarity in polarities:
        for particle in particles:
            outfile =  datadir+'MISID/SameSign/'+particle+'_sample_'+polarity+suffix[sample]
            print('>>>>>   Processing:          ', particle, ' ', polarity )
            print('Sample = ', sample)
            print('Output file name = ', outfile)
            AddMISIDweightNoCF(files[polarity], outfile, sample,polarity, particle)




