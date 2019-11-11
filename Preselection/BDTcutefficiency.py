import ROOT as r
import numpy as np
from array import array

from ROOT import gSystem
gSystem.Load('libRooFit')
from ROOT import RooStats
from ROOT import RooFit as RF
from ROOT import RooRealVar, RooGaussian, RooDataSet, RooArgList, RooTreeData, RooDataHist

from AddSweights import *
#from rootpy.interactive import wait

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
files = {'Data':datadir + 'Data/Lb_Data_MagUp_reduced.root','MC':datadir+'MC/Lb_Lctaunu_PID_reduced.root'}

def Fit(h,Lc_M):
    h_LcM = RooDataHist('h_LcM','h_LcM',RooArgList(Lc_M),h)

    #Create a workspace named 'w' with the model to fit the Lc_Mass peak
    nsig = RooRealVar("nsig","N signal evts",100000,0,5000000)
    nbkg = RooRealVar("nbkg","N bkg evts",400000,0,1E10)

    w = CreateMassFitModel(Lc_M,nsig,nbkg)
    model = w.pdf("model")
    model.fitTo(h_LcM,RF.Extended(True))
    return model, w, h_LcM

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

def GetInitialSignalBkg(fname):
    f = r.TFile(fname,'read')
    t = f.Get('DecayTree')

    LcM_range = [2230,2330]
    Lc_M = RooRealVar('Lc_M','#Lambda_{c} mass',LcM_range[0],LcM_range[1])
    h = r.TH1F('h','#Lambda_{c} mass peak',500,2230,2330)
    t.Draw('Lc_M>>h')
    h = r.gPad.GetPrimitive('h')
    h = ScaleHisto(h,1)

    model, w, h_LcM = Fit(h,Lc_M)
    nsig = w.var('nsig')
    nbkg = w.var('nbkg')

    c = r.TCanvas()
    frame = Lc_M.frame(RF.Title('#Lambda_{c} mass peak distribution'))
    h_LcM.plotOn(frame)
    model.plotOn(frame)
    model.paramOn(frame, RF.Layout(0.6,0.9,0.85))
    model.plotOn(frame, RF.Components('bkg'),RF.LineColor(r.kRed))
    frame.Draw()
    #c.Print()
    return nsig.getValV(), nbkg.getValV()


def EfficiencySignal(fname, bdtcut):
    f = r.TFile(fname,'read')
    t = f.Get('DecayTree')
    
    LcM_range = [2230,2330]
    Lc_M = RooRealVar('Lc_M','#Lambda_{c} mass',LcM_range[0],LcM_range[1])

    h0 = r.TH1F('h0','#Lambda_{c} mass peak',200,2230,2330)
    t.Draw('Lc_M>>h0','Lc_BKGCAT<30.')
    h0 = r.gPad.GetPrimitive('h0')


    model, w, h_LcM = Fit(h0,Lc_M)
    nsig_nocut = w.var('nsig')
    nbkg_nocut = w.var('nbkg')


    h1 = r.TH1F('h1','#Lambda_{c} mass peak with BDT cut',200,2230,2330)
    t.Draw('Lc_M>>h1','Lc_BKGCAT<30.&&bdt>'+str(bdtcut))
    h1 = r.gPad.GetPrimitive('h1')


    model1, w1, h_LcM1= Fit(h1,Lc_M)
    nsig_cut = w1.var('nsig')
    nbkg_cut = w1.var('nbkg')

    return nsig_cut.getValV()/nsig_nocut.getValV()


    '''
    c = r.TCanvas('c','c')
    frame = Lc_M.frame(RF.Title('#Lambda_{c} mass peak distribution'))
    h_LcM.plotOn(frame)
    model.plotOn(frame)
    model.paramOn(frame, RF.Layout(0.6,0.9,0.85))
    model.plotOn(frame, RF.Components('bkg'),RF.LineColor(r.kRed))
    frame.Draw()
    c.Print()
    wait()
    '''

def EfficiencyBkg(fname, bdtcut):
    f = r.TFile(fname,'read')
    t = f.Get('DecayTree')
    
    entries_nocut = t.Draw("Lc_M","Lc_M<2260. || Lc_M>2310.")
    entries_cut = t.Draw("Lc_M","(Lc_M<2260. || Lc_M>2310.)&&bdt>"+str(bdtcut))

    return entries_cut*1./entries_nocut





if __name__ == '__main__':
    nsig, nbkg = GetInitialSignalBkg(files['Data'])
    eff_s, eff_b, signif, signif_1, bdtcut = array( 'd' ), array( 'd' ), array( 'd' ),array( 'd' ),array( 'd' )
 

    n=0
    for i in np.arange(0,1,0.05):
        n=n+1
        bdtcut.append(i)
        print('bdtcut = ',i)
        eff_s.append(EfficiencySignal(files['MC'], i))
        eff_b.append(EfficiencyBkg(files['Data'],i))
        signif.append(eff_s[n-1]*nsig/r.TMath.Sqrt(nsig*eff_s[n-1]+nbkg*eff_b[n-1]))
        signif.append(eff_s[n-1]*nsig/r.TMath.Sqrt(nsig*eff_s[n-1]+nbkg*eff_b[n-1]))

        #print(eff_s, eff_b) 
    c = r.TCanvas('c','c',500,500)
    gr = r.TGraph( n, bdtcut, eff_s )
    gr1 = r.TGraph( n, bdtcut, eff_b )
    gr2 = r.TGraph( n, bdtcut, signif )
    gr.SetLineColor( 2 )
    gr.SetLineWidth( 4 )
    gr.SetMarkerColor( 2 )
    gr.SetMarkerStyle( 21 )
    gr.GetYaxis().SetRangeUser(0,1)
    gr.GetXaxis().SetTitle( 'BDT_Cut' )
    gr.GetYaxis().SetTitle( 'Eff_s' )
    gr1.SetLineColor( 4 )
    gr1.SetLineWidth( 4 )
    gr1.SetMarkerColor( 4 )
    gr1.SetMarkerStyle( 21 )
    gr1.GetXaxis().SetTitle( 'BDT_Cut' )
    gr1.GetYaxis().SetTitle( 'Eff_bkg' )
    gr2.SetLineColor( 1 )
    gr2.SetLineWidth( 4 )
    gr2.SetMarkerColor( 1 )
    gr2.SetMarkerStyle( 21 )
    gr2.GetXaxis().SetTitle( 'BDT_Cut' )
    gr2.GetYaxis().SetTitle( 'S/#sqrt(S+B)' )
    text = r.TPaveText(.7,.87,.8,.92);
    text.AddText("#epsilon_{s}") 
    text.GetListOfLines().Last().SetTextColor(r.kRed)
    text1 = r.TPaveText(.7,.15,.8,.2);
    text1.AddText("#epsilon_{b}") 
    text1.GetListOfLines().Last().SetTextColor(r.kBlue)
    text2 = r.TPaveText(.7,.45,.8,.5);
    text2.AddText("S/#Sqrt(S+B)") 
    text2.GetListOfLines().Last().SetTextColor(r.kBlack)


    gr.Draw( 'ACP' )
    gr1.Draw('CP')
    gr2.Draw('CP')
    '''
    text.Draw('same')
    text1.Draw('same')
    text2.Draw('same')
'''

