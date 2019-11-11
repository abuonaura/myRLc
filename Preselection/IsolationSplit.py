import ROOT as r
import os,sys,getopt,time
import numpy as np
import argparse
from AddSweights import *

parser = argparse.ArgumentParser(description='Split samples accoriding to isolation requests + add sweights')
parser.add_argument('datatype',choices=['Data','FakeMu','FakeMuSS','DataSS','all'], help = 'which data sample we want to run on')
parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on')
parser.add_argument('category', choices=['full','iso','all'], help ='which isolation category we want to process')

args = parser.parse_args()

datatype =args.datatype
category = args.category
polarity=args.polarity

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
polarities = ['MagUp','MagDown']

ISOBDTcut = 0.35
ISOBDT2cut = 0.2
isocut = {'full': '', 'iso':"Lb_ISOLATION_BDT<"+str(ISOBDTcut)} 

suffix = {'full':'.root', 'iso':'_iso.root'}

#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)



def EvaluateCutEfficiency(t,ot):
    start = t.GetEntries()
    end = ot.GetEntries()
    eff = end*1./start
    return eff


def SplitSample(fname,cat):
    f = r.TFile(fname,'read')
    t = f.Get('DecayTree')
    tl = f.Get('LumiTuple')

    ofname = fname[0:-5]+suffix[cat]
    
    of = r.TFile(ofname,'recreate')
    ot = r.TTree('DecayTree','DecayTree')
    otl =r.TTree('LumiTuple','LumiTuple')

    print(' >>> Copying tree')
    print(' Applied cuts: '+isocut[cat])
    ot = t.CopyTree(isocut[cat])
    otl = tl.CloneTree()

    ot.Write()
    
    eff = EvaluateCutEfficiency(t,ot)
    print('The cut selects %f of the original events' %(eff))

    otl.Write()
    of.Close()
    f.Close()
    return ofname


if __name__== "__main__":
    if datatype=='all':
        datatypes = ['Data','FakeMu','FakeMuSS','DataSS']
    else:
        datatypes=[datatype]
    if category=='all':
        categories = ['full','iso','Kenriched']
    else:
        categories = [category]
    if polarity!='all':
        polarities=[polarity]

    for dtype in datatypes:
        print('>>>>   Processing %s sample' %(dtype))
        for pol in polarities:
            fname = datadir+'Data/Lb_'+dtype+'_'+pol+'_reduced_preselected.root'
            print('- Polarity %s' %pol)
            print('- Input file: %s' %fname)
            for cat in categories:
                print('-category: %s' %cat)
                ofname = SplitSample(fname, cat)
                print('- created file: %s',ofname)
                swfname =ComputeSweights(ofname, dtype, pol,cat)
                tname = 'DecayTree'
                variable = 'Lc_M'
                splotVariable(swfname, tname, dtype, pol,cat,variable, 50, 2230, 2330, 'm_{#Lambda_{c}}')





