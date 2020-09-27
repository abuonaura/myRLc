'''
author: A. Buonaura
Date: 10th June 2020
Comment: This scripts compare the ISOLATION BDT distributions for Lb_LcDs and Lb_Lctaunu and try to define a cut to separate signal and bkg region.
'''

import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas, TMath, TLegend
from array import array

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'

datafiles={'Signal':{'MagUp': filedir+''}}
