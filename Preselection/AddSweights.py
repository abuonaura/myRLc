import ROOT as r
import numpy as np
import os

from ROOT import gSystem
gSystem.Load('libRooFit')
from ROOT import RooStats
from ROOT import RooFit as RF
from ROOT import RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData

suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

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

def splotVariable(fname, tname, datatype, polarity,sample,variable, Nbin, min, max, Xaxis='', Yaxis=''):
    f = r.TFile(fname,"READ")
    #if os.path.isfile(fname):
    print('File found')
    t = f.Get(tname)
        
    var = np.zeros(1, dtype=float)
    sw = np.zeros(1, dtype=float)
        
    histname = variable + '_sweighted'
    h = r.TH1F(histname,histname,Nbin,min,max)
    h.GetXaxis().SetTitle(Xaxis)
    h.Sumw2()
        
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        h.Fill(getattr(t,variable),t.sw_sig)
        h.SetDirectory(0)

    c = r.TCanvas()
    h.Draw()
    suffix_img = suffix[sample][0:-5]+'.png'
    c.SaveAs('plots/'+histname+'_'+datatype+'_'+polarity+suffix_img)
    return 

def PlotMassFit(model, ds, Lc_M,datatype,polarity,sample):
    c = r.TCanvas()
    frame = Lc_M.frame(RF.Title('#Lambda_{c} mass peak distribution'))
    ds.plotOn(frame)
    model.plotOn(frame)
    model.paramOn(frame, RF.Layout(0.6,0.9,0.85))
    model.plotOn(frame, RF.Components('bkg'),RF.LineColor(r.kRed))
    frame.Draw()
    c.Draw()
    suffix_img = suffix[sample][0:-5]+'.png'
    c.SaveAs('FitResults/Fit_'+datatype+'_'+polarity+suffix_img)
    return 

def ComputeSweights(fname, datatype, polarity,sample):
    f = r.TFile(fname,'update')
    t = f.Get('DecayTree')
    ofname = fname[0:-5]+'_sw.root'
    #of = r.TFile(ofname,'recreate')
    #ot = r.TTree('DecayTree','DecayTree')
    
    #print(' ... Copying tree to add sweights ...')
    #ot = t.CloneTree(0)

    LcM_range = [2230,2330]
    Lc_M = r.RooRealVar('Lc_M','#Lambda_{c} mass',LcM_range[0],LcM_range[1])
    ds = GetDataset(t,Lc_M)

    #Create a workspace named 'w' with the model to fit the Lc_Mass peak
    nsig = r.RooRealVar("nsig","N signal evts",100000,0,1000000)
    nbkg = r.RooRealVar("nbkg","N bkg evts",400000,0,1000000)

    w = CreateMassFitModel(Lc_M,nsig,nbkg)
    model = w.pdf("model")
    result = model.fitTo(ds,RF.Extended(True),RF.Save())


    nsig = w.var('nsig')
    nbkg = w.var('nbkg')

    sData =r.RooStats.SPlot("sData","An SPlot",ds, model, r.RooArgList(nsig, nbkg))

    ds_w = GetWeightedDS(ds)

    Lc_M_data = np.zeros(1, dtype=float)
    sw_bkg = np.zeros(1, dtype=float)
    sw_sig = np.zeros(1, dtype=float)
    
    if t.FindBranch('sw_bkg') and t.FindBranch('sw_sig'):
        sw_bkg_b = t.GetBranch('sw_bkg')
        sw_sig_b = t.GetBranch('sw_sig')
        sw_bkg_b.Reset()
        sw_sig_b.Reset()
    sw_bkg_b = t.Branch("sw_bkg",sw_bkg,"sw_bkg/D")
    sw_sig_b = t.Branch("sw_sig",sw_sig,"sw_sig/D")

    #print('<<<<< Created file ', ofname)
    
    for i in range(ds_w.numEntries()):
        t.GetEntry(i)
        sw_bkg[0]= ds_w.get(i).getRealValue("nbkg_sw")
        sw_sig[0] = ds_w.get(i).getRealValue("nsig_sw")
        sw_bkg_b.Fill()
        sw_sig_b.Fill()

    t.Write()
    f.Close()
    
    os.system('cp '+fname+ ' '+ofname)
    #Save fit results to file
    of_fit = r.TFile('FitResults/FitResults_'+datatype+'_'+polarity+suffix[sample],'Recreate')
    result.Write()
    of_fit.Close()
    
    PlotMassFit(model, ds, Lc_M,datatype,polarity,sample)
    return ofname



