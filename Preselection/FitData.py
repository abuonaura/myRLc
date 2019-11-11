import ROOT as r
import numpy as np

from ROOT import gSystem
gSystem.Load('libRooFit')
from ROOT import RooStats
from ROOT import RooFit as RF
from ROOT import RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData
from ROOT import TFile, TCanvas

import os,sys,getopt,time
import argparse

parser = argparse.ArgumentParser(description='Evaluate signal yield')
parser.add_argument('datatype',choices=['Data','FakeMu','FakeMuSS','DataSS','all'], help = 'which data sample we want to run on')
parser.add_argument('category', choices=['full','iso','Kenriched','all'], help ='which isolation category we want to process')

args = parser.parse_args()

datatype =args.datatype
category = args.category

suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}
polarities = ['MagUp','MagDown']

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/Data/'
filename = {'MagUp':datadir+'Lb_Data_MagUp_reduced_preselected.root','MagDown':datadir+'Lb_Data_MagDown_reduced_preselected.root'}

if __name__== "__main__":
    if datatype=='all':
        datatypes = ['Data','FakeMu','FakeMuSS','DataSS']
    else:
        datatypes=[datatype]
    if category=='all':
        categories = ['full','iso','Kenriched']
    else:
        categories = [category]

    for dtype in datatypes:
        print('>>>>   Processing %s sample' %(dtype))
        for polarity in polarities:
            print('- Polarity %s' %polarity)
            for cat in categories:
                print('-category: %s' %cat)
                fname = datadir+'Lb_'+dtype+'_'+polarity+'_reduced_preselected'+suffix[cat]
                print('- Input file: %s' %fname)
                f = TFile(fname,'READ')
                t = f.Get('DecayTree')

                LcM_range = [2230,2330]
                Lc_M = r.RooRealVar('Lc_M','#Lambda_{c} mass',LcM_range[0],LcM_range[1])

                ds = r.RooDataSet("ds","ds",t,r.RooArgSet(Lc_M))

                #Create a workspace named 'w' with the model to fit the Lc_Mass peak
                if cat=='Kenriched':
                    nsig = r.RooRealVar("nsig","N signal evts",20000,0,30000)
                    nbkg = r.RooRealVar("nbkg","N bkg evts",10000,0,30000)
                else:
                    nsig = r.RooRealVar("nsig","N signal evts",500000,0,1000000)
                    nbkg = r.RooRealVar("nbkg","N bkg evts",400000,0,500000)


                w = r.RooWorkspace()

                #Create two Gaussian PDFs g1(x,mean1,sigma) anf g2(x,mean2,sigma) and their parameters
                mean = r.RooRealVar("mean","mean of gaussians",2290,2270,2310) 
                if cat=='Kenriched':
                    sigma1 = r.RooRealVar("sigma1","width of gaussians",2,0,7)
                else:
                    sigma1 = r.RooRealVar("sigma1","width of gaussians",3,0,10) 
                if cat=='Kenriched':
                    sigma2 = r.RooRealVar("sigma2","width of gaussians",7,0,15) 
                else:
                    sigma2 = r.RooRealVar("sigma2","width of gaussians",15,0,30) 
                sig1 = r.RooGaussian("sig1","Signal component 1",Lc_M,mean,sigma1)
                sig2 = r.RooGaussian("sig2","Signal component 2",Lc_M,mean,sigma2)

                #Build exponential PDF
                l = r.RooRealVar("l", "slope of expo",-0.1, -1., 0.);
                bkg = r.RooExponential("bkg", " bkg with exponential PDF",Lc_M,l)

                #Sum the signal components into a composite signal p.d.f.
                sig1frac = r.RooRealVar("sig1frac","fraction of component 1 in signal",0.5,0.,1.)
                sig = r.RooAddPdf("sig","Signal",sig1,sig2,sig1frac)

                yields = r.RooArgList(nsig,nbkg)
                functions = r.RooArgList(sig,bkg)
                #Sum the composite signal and background
                model= r.RooAddPdf("model","g1+g2+a",functions,yields)

                getattr(w,'import')(model)

                model = w.pdf("model")

                result = model.fitTo(ds,RF.Extended(True),RF.Save())

                #Save fit results to file
                of_fit = TFile('FitResults_Data_'+polarity+'.root','Recreate') 
                result.Write()
                of_fit.Close()


                #Plot Fit results
                c = r.TCanvas()
                frame = Lc_M.frame(RF.Title('#Lambda_{c} mass peak distribution'))
                ds.plotOn(frame)
                model.plotOn(frame)
                model.plotOn(frame, RF.Components('bkg'),RF.LineColor(r.kRed))
                model.paramOn(frame, RF.Layout(0.6,0.9,0.85))
                frame.Draw()
                c.Draw()
                c.SaveAs('plots/NewFit_'+dtype+'_'+polarity+'_'+cat+'.png')
