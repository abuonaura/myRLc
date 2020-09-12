'''
Author: A. Buonaura
Date: 10th July 2019

'''

import ROOT as r
import numpy as np
import os, sys, getopt, time

LbMass=5620.2
datadir = '/disk/lhcb_data2/RLcMuonic2016/'
polarities = ['MagUp','MagDown']

#Lb_DataSS_MagDown_sw_noMISID.root

suffix = {'iso':'_iso.root','Kenriched':'_Kenr.root','Lcpipi':'_Lcpipi.root'}

if __name__ == '__main__':
    reweight=0
    opts, args = getopt.getopt(sys.argv[1:], "",["iso","Kenriched","Lcpipi","reweight"])
    print (opts,args)
    for o, a in opts:
        if o in ("--iso",):
            sample = 'iso'
        if o in ("--Kenriched",):
            sample = 'Kenriched'
        if o in ("--Lcpipi",):
            sample = 'Lcpipi'
        if o in ("--reweight",):
            reweight=1

    for polarity in polarities:
        fname = datadir + 'Data/Lb_DataSS_'+polarity+'_sw_noMISID'+suffix[sample]
        ofname = datadir + 'CombinatorialBkg/CombinatorialBkg_'+polarity+suffix[sample]
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')

        of = r.TFile(ofname,'recreate')
        ot = r.TTree('DecayTree','DecayTree')
        ot = t.CloneTree(0)

        w_comb = np.zeros(1,dtype=float) 
        ot.Branch("w_comb",w_comb,"w_comb/D")
        fw = r.TF1('fw','[0]*exp([1]*x+[2])',5000,10000)
        fw.FixParameter(0,8.62114e-01)
        fw.FixParameter(1,-1.36285e-04)
        fw.FixParameter(2,1.09663e+00)


        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_M<LbMass:
                w_comb[0] = fw.Eval(t.Lb_M)
                ot.Fill()

        of.Write()
        of.Close()

        f.Close()

    h_MagUp = r.TH3F('h_MagUp','; q2^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)

    h_MagDown = r.TH3F('h_MagDown','; q2^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)

    fname_Up = datadir + 'CombinatorialBkg/CombinatorialBkg_MagUp'+suffix[sample]
    fUp = r.TFile(fname_Up,'READ')
    tUp = fUp.Get('DecayTree')
    for i in range(tUp.GetEntries()):
        tUp.GetEntry(i)
        if reweight==1:
            h_MagUp.Fill(tUp.FitVar_q2_mLc/1E6,tUp.FitVar_El_mLc,tUp.FitVar_Mmiss2_mLc/1E6,tUp.sw_sig*tUp.w_MISID*tUp.w_comb)
        else:
            h_MagUp.Fill(tUp.FitVar_q2_mLc/1E6,tUp.FitVar_El_mLc,tUp.FitVar_Mmiss2_mLc/1E6,tUp.sw_sig*tUp.w_MISID)

    fname_Down = datadir + 'CombinatorialBkg/CombinatorialBkg_MagDown'+suffix[sample]
    fDown = r.TFile(fname_Down,'READ')
    tDown = fDown.Get('DecayTree')
    for i in range(tDown.GetEntries()):
        tDown.GetEntry(i)
        if reweight==1:
            h_MagDown.Fill(tDown.FitVar_q2_mLc/1E6,tDown.FitVar_El_mLc,tDown.FitVar_Mmiss2_mLc/1E6,tDown.sw_sig*tDown.w_MISID*tDown.w_comb)
        else:
            h_MagDown.Fill(tDown.FitVar_q2_mLc/1E6,tDown.FitVar_El_mLc,tDown.FitVar_Mmiss2_mLc/1E6,tDown.sw_sig*tDown.w_MISID)


    h_CombBkg = r.TH3F('h_CombBkg','; q2^{2} (Gev^{2}); E_{l} (MeV^{2}); Mmiss^{2} (GeV^{2})',4,-2,14,10,0,2600,10,-2,14)
    h_CombBkg.Add(h_MagUp, h_MagDown)

    c = r.TCanvas('c','',1500,500)
    c.Divide(3,1)
    c.cd(1)
    h_CombBkg.ProjectionX().Draw()
    c.cd(2)
    h_CombBkg.ProjectionY().Draw()
    c.cd(3)
    h_CombBkg.ProjectionZ().Draw()

    print('Combinatorial Background Yield: ', h_CombBkg.Integral())
