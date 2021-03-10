#!/bin/python

from ROOT import TFile, TTree, TH2D, TH1, kFALSE
import uproot4
import uproot3 as uproot3
import numpy as np
from store_q2costhmu import *


def storeeff2D(row, histg, limits, q2factor):
    q2val  = row[0]*q2factor
    cthlval= row[1]
    glbbin = histg.FindBin(q2val, cthlval) #global bin number
    eff    = histg.GetBinContent(glbbin)
    xmin, xmax, ymin, ymax = limits
    cond   = (q2val < xmin or q2val > xmax or cthlval < ymin or cthlval > ymax)
    if cond: eff = 0.
    return eff

def AddFFcorr(infname, intreename, outfname, outtreename, Lcstate, leptname, q2True_branchname, costhlTrue_branchname, nentries_to_read = 1000000000, chunksize = 10000):

    TH1.AddDirectory(kFALSE)

    perfname = None
    q2factor = None
    if Lcstate == 'Lc':
        perfname = './CorrectionTables/LcFFratios.root'
        q2factor = 1.
    elif (Lcstate == 'Lc2595' or Lcstate == 'Lc2625'):
        perfname = './CorrectionTables/LcstFFratios.root'
        q2factor = 1e-6
    else:
        raise Exception('Lc state not recognised', Lcstate)
    
    if leptname != 'mu' and leptname != 'tau':
        raise Exception('Lepton name not recognised', leptname)

    print('Using the histname', Lcstate+leptname+"_ratio")

    #variables to get from file
    varsdf    = ['runNumber', 'eventNumber'] 
    varsdf   += ['Lb_TRUEP_X','Lb_TRUEP_Y','Lb_TRUEP_Z','Lb_TRUEP_E']
    varsdf   += ['Lc_TRUEP_X','Lc_TRUEP_Y','Lc_TRUEP_Z','Lc_TRUEP_E']
    varsdf   += ['Lb_True'+leptname.capitalize()+'_PX', 'Lb_True'+leptname.capitalize()+'_PY', 'Lb_True'+leptname.capitalize()+'_PZ', 'Lb_True'+leptname.capitalize()+'_PE']
    varsdf   += ['Lb_TrueNeutrino_PX', 'Lb_TrueNeutrino_PY', 'Lb_TrueNeutrino_PZ', 'Lb_TrueNeutrino_PE']

    File      = TFile.Open(perfname, "read")
    Histg     = File.Get(Lcstate+leptname+"_ratio")
    perfHist  = Histg.Clone(Lcstate+leptname+"_rationew")
    File.Close()
    Xmin  = perfHist.GetXaxis().GetXmin(); Xmax  = perfHist.GetXaxis().GetXmax()
    Ymin  = perfHist.GetYaxis().GetXmin(); Ymax  = perfHist.GetYaxis().GetXmax()
    Limits= (Xmin, Xmax, Ymin, Ymax)
    print(Limits, perfHist.Integral())

    #variables to store in the new ttree
    varstoStore = {'runNumber'    : np.int, 
                   'eventNumber'  : np.int, 
                   'Event_FFcorr': np.float64,
                   costhlTrue_branchname: np.float64,
                   q2True_branchname : np.float64}

    aliases = {}
    #create a new rootfile
    with uproot3.recreate(outfname) as f:
        f[outtreename] = uproot3.newtree(varstoStore)

        #loop over the old rootfile chunkwise
        events_read = 0
        if chunksize >= nentries_to_read: chunksize = nentries_to_read
        for df_data in uproot4.iterate(infname+':'+intreename,varsdf,aliases=aliases,cut=None,library="pd",step_size=chunksize): 
            if events_read >= nentries_to_read: break

            #Compute q2 and cosThetaL
            pxl  = df_data['Lb_True'+leptname.capitalize()+'_PX']; pxnu = df_data['Lb_TrueNeutrino_PX']
            pyl  = df_data['Lb_True'+leptname.capitalize()+'_PY']; pynu = df_data['Lb_TrueNeutrino_PY']
            pzl  = df_data['Lb_True'+leptname.capitalize()+'_PZ']; pznu = df_data['Lb_TrueNeutrino_PZ']
            pel  = df_data['Lb_True'+leptname.capitalize()+'_PE']; penu = df_data['Lb_TrueNeutrino_PE']
            if (Lcstate == 'Lc2595' or Lcstate == 'Lc2625'):
                #this should be Lcstar momentum
                pxlc = df_data['Lb_TRUEP_X'] - pxl - pxnu
                pylc = df_data['Lb_TRUEP_X'] - pyl - pynu
                pzlc = df_data['Lb_TRUEP_X'] - pzl - pznu
                pelc = df_data['Lb_TRUEP_X'] - pel - penu
            elif Lcstate == 'Lc':
                pxlc = df_data['Lc_TRUEP_X']
                pylc = df_data['Lc_TRUEP_Y']
                pzlc = df_data['Lc_TRUEP_Z']
                pelc = df_data['Lc_TRUEP_E']

            PLc_lab   = LorentzVector(Vector(pxlc, pylc, pzlc), pelc) #Format of LorentzVector(Vector(X,Y,Z), E)
            Pl_lab    = LorentzVector(Vector(pxl , pyl , pzl ), pel )
            PNu_lab   = LorentzVector(Vector(pxnu, pynu, pznu), penu)
            PLb_lab   = PLc_lab + Pl_lab + PNu_lab
            qsq, cthl = return_phasespace(PLb_lab, PLc_lab, Pl_lab)
            #print(qsq,cthl)
            df_data[q2True_branchname] = qsq
            df_data[costhlTrue_branchname] = cthl

            #get the corrections
            applyvars = [q2True_branchname, costhlTrue_branchname] #has to be in correct order like in histogram
            df_data['Event_FFcorr'] = df_data[applyvars].apply(storeeff2D, args=[perfHist, Limits, q2factor], axis=1)

            #get only the things that need to be stored and write them to the file
            branch_dict = {vartostore: df_data[vartostore].to_numpy() for vartostore in list(varstoStore.keys())}
            f[outtreename].extend(branch_dict)
            events_read += df_data.shape[0]
            print('Events read', events_read)
