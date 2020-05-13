import ROOT as r
import math as m
import sys, os

from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.externals import joblib
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
from array import array
rando = r.TRandom3()

def random_smearing(P,PT,realET,mom_high,mom_low,mombins,pt_high,pt_low,ptbins,momPThists):
# choosing the P-Pt bin/histo for this particle
  Pbin = int(m.floor(P/((mom_high-mom_low)/mombins)))         # m.floor giving back the biggest round number
  PTbin = int(m.floor(PT/((pt_high-pt_low)/ptbins)))
  if Pbin < 0: Pbin=0
  if PTbin < 0: PTbin=0
  if Pbin >= mombins: Pbin=mombins-1
  if PTbin >= ptbins: PTbin=ptbins-1
# taking the right momentum region histo of calorimeter uncertainties 
  momPThist = momPThists[Pbin][PTbin]
# take a random value for the uncertainties from the histo
  rando = momPThist.GetRandom()
# realET is the energy value from the tracker
  smearedET = realET*(1-rando)
  if smearedET < 0:
    smearedET=0
  if smearedET > 6100:  # it is limitation from HCAL. It will be overloaded wth bigger values
    smearedET=6100 
  return smearedET 


def missingFraction(rdiff,region,clusters):
  hmissing = clusters['missing']['outer']
  hshared = clusters['shared']['outer']
  if region == 1:
    hmissing = clusters['missing']['inner']
    hshared = clusters['shared']['inner']

  fracmissing = hmissing.GetBinContent(hshared.FindBin(rdiff))/hmissing.GetBinContent(10)
  #fracmissing=1.0

  #print(rdiff,region,hmissing.GetBinContent(hshared.FindBin(rdiff)),hmissing.GetBinContent(10))
  return fracmissing

def isShared(rdiff,region,clusters):
  hshared = clusters['shared']['outer']
  if region == 1:
    hshared = clusters['shared']['inner']
  fracshared = hshared.GetBinContent(hshared.FindBin(rdiff))
  #print('frac shared',fracshared,'rdiff',rdiff)
  #rando = r.TRandom3()
  if rando.Uniform() < fracshared:
    return True
  else:
    return False
def main(MC_fileName = "../../../fromSimone/L0HLT1emulationVars/Lb_Lcmunu_MagUp.root"):
    
    mombins = 10;
    ptbins = 6;
    pt_low = 0;
    pt_high = 15000;
    mom_low = 0;
    mom_high = 1e5;

    TOSthreshold = 3744 #3500 ; changed threshhold because of the new MC version

    momhistsfile = r.TFile("momhists.root")

    momPThists = []

    for i in range(mombins):
        PThists = []
        for j in range(ptbins):
            hname= "hdiff_"+str(i)+"_"+str(j)
            PThists.append(momhistsfile.Get(hname)) 
        momPThists.append(PThists)


    rando = r.TRandom3()
    fcluster_hists = r.TFile.Open("two_particle_HCAL_cluster_hists.root")
    clusters={'shared':{'outer':0,'inner':0},'missing':{'outer':0,'inner':0}}
    for iC in clusters.keys():
        for iW in clusters[iC].keys():
            clusters[iC][iW] = fcluster_hists.Get(iC+"_with_radial_"+iW)




    fin = r.TFile.Open(MC_fileName)
    tin = fin.Get('tupleout/DecayTree')
    partTypes = ['K','p','pi']#,'Lc']
    physObs = ['_L0Calo_HCAL_xProjection','_L0Calo_HCAL_yProjection','_L0Calo_HCAL_realET','_L0Calo_HCAL_region','_P','_PT']
    nentries = tin.GetEntries()
    tin.SetBranchStatus('*',0)

    for iP in list(itertools.product(partTypes,physObs)):
        if tin.FindBranch(iP[0]+iP[1]):
            tin.SetBranchStatus(iP[0]+iP[1],1)
    tin.SetBranchStatus("Lc_PT",1)
    tin.SetBranchStatus("Lc_P",1)
    fout = r.TFile.Open(MC_fileName[:-5]+'_wL0TOSEmulation.root','RECREATE')
    tout = r.TTree("DecayTree", "DecayTree")

    inputarray = np.zeros((nentries,8))
    Lc_ET_emulated_array = np.zeros(nentries)

    print("Processing events ...")
    for i,ev in enumerate(tin):
        #print(i)
        if i%100 == 0:
                print(i)
                #if i==100: break
        daughtParams = {}
        daughtCalc = {}
        for iP in list(itertools.product(partTypes,physObs)):
            if tin.FindBranch(iP[0]+iP[1]):
                daughtParams[iP[0]+iP[1]] = ev.GetLeaf(iP[0]+iP[1]).GetValue()
                #print(iP[0]+iP[1])
                
        for iP in range(3):
            daughtCalc[partTypes[iP]+'_smeared'] = random_smearing(daughtParams[partTypes[iP]+'_P'],
                                                                    daughtParams[partTypes[iP]+'_PT'],
                                                                    daughtParams[partTypes[iP]+'_L0Calo_HCAL_realET'],mom_high,mom_low,mombins,pt_high,pt_low,ptbins,momPThists
                                                                    )
            if daughtParams[partTypes[iP]+'_L0Calo_HCAL_region'] == -1:
                daughtCalc[partTypes[iP]+'_smeared'] = 0
            daughtCalc[partTypes[iP]+'_ET'] = daughtCalc[partTypes[iP]+'_smeared']

        for iP in range(3):
            for iAP in range(1,3):
                if iP<iAP:
                    daughtCalc['rdiff_'+partTypes[iP]+partTypes[iAP]] = m.sqrt( 
                        (  daughtParams[partTypes[iP]+'_L0Calo_HCAL_xProjection'] - daughtParams[partTypes[iAP]+'_L0Calo_HCAL_xProjection']
                        )**2
                    +  
                        (  daughtParams[partTypes[iP]+'_L0Calo_HCAL_yProjection'] - daughtParams[partTypes[iAP]+'_L0Calo_HCAL_yProjection']
                        )**2
                    )
                        
        #             if daughtParams[partTypes[iP]+'_L0Calo_HCAL_region'] == daughtParams[partTypes[iAP]+'_L0Calo_HCAL_region']:

        #                 if isShared(daughtCalc['rdiff_'+partTypes[iP]+partTypes[iAP]],daughtParams[partTypes[iP]+'_L0Calo_HCAL_region'],clusters):
        #                         daughtCalc[partTypes[iP]+'_ET'] = daughtCalc[partTypes[iP]+'_smeared'] + daughtCalc[partTypes[iAP]+'_smeared']
        #                         daughtCalc[partTypes[iAP]+'_ET'] = daughtCalc[partTypes[iP]+'_ET']
        #                         #print('do missing corr',i)
        #                         daughtCalc[partTypes[iP]+'_ET'] = daughtCalc[partTypes[iP]+'_ET']*missingFraction( daughtCalc['rdiff_'+partTypes[iP]+partTypes[iAP]], daughtParams[partTypes[iP]+'_L0Calo_HCAL_region'],clusters )
        #                         daughtCalc[partTypes[iAP]+'_ET'] = daughtCalc[partTypes[iAP]+'_ET']*missingFraction( daughtCalc['rdiff_'+partTypes[iP]+partTypes[iAP]],
        #                                                                             daughtParams[partTypes[iP]+'_L0Calo_HCAL_region'],clusters)
        # daughtCalc[partTypes[0]+'_ET'] = daughtCalc[partTypes[0]+'_smeared'] 
        # daughtCalc[partTypes[1]+'_ET'] = daughtCalc[partTypes[1]+'_smeared']
        # daughtCalc[partTypes[2]+'_ET'] = daughtCalc[partTypes[2]+'_smeared']
        # if (i == 1):
        #   print(daughtCalc['rdiff_Kminuspiplus'],daughtCalc['rdiff_pipluspiplus0'],daughtCalc['rdiff_Kminuspiplus0'])

        if daughtParams[partTypes[0]+'_L0Calo_HCAL_region'] == daughtParams[partTypes[1]+'_L0Calo_HCAL_region']:
          if isShared(daughtCalc['rdiff_'+partTypes[0]+partTypes[1]],daughtParams[partTypes[0]+'_L0Calo_HCAL_region'],clusters):
            # print('before')
            # print(daughtCalc[partTypes[0]+'_ET'],daughtCalc[partTypes[1]+'_ET'])
            daughtCalc[partTypes[0]+'_ET'] = (daughtCalc[partTypes[0]+'_smeared']+daughtCalc[partTypes[1]+'_smeared'])*missingFraction( daughtCalc['rdiff_'+partTypes[0]+partTypes[1]], daughtParams[partTypes[0]+'_L0Calo_HCAL_region'],clusters )
            daughtCalc[partTypes[1]+'_ET'] = daughtCalc[partTypes[0]+'_ET']
            # print('after')
            # print(daughtCalc[partTypes[0]+'_ET'],daughtCalc[partTypes[1]+'_ET'])

        if daughtParams[partTypes[1]+'_L0Calo_HCAL_region'] == daughtParams[partTypes[2]+'_L0Calo_HCAL_region']:
          if isShared(daughtCalc['rdiff_'+partTypes[1]+partTypes[2]],daughtParams[partTypes[1]+'_L0Calo_HCAL_region'],clusters):
            daughtCalc[partTypes[1]+'_ET']+=daughtCalc[partTypes[2]+'_ET']*missingFraction( daughtCalc['rdiff_'+partTypes[1]+partTypes[2]], daughtParams[partTypes[1]+'_L0Calo_HCAL_region'],clusters )
            daughtCalc[partTypes[2]+'_ET'] = (daughtCalc[partTypes[1]+'_smeared']+daughtCalc[partTypes[2]+'_smeared'])*missingFraction( daughtCalc['rdiff_'+partTypes[1]+partTypes[2]], daughtParams[partTypes[1]+'_L0Calo_HCAL_region'],clusters )

        if daughtParams[partTypes[0]+'_L0Calo_HCAL_region'] == daughtParams[partTypes[2]+'_L0Calo_HCAL_region']:
          if isShared(daughtCalc['rdiff_'+partTypes[0]+partTypes[2]],daughtParams[partTypes[0]+'_L0Calo_HCAL_region'],clusters):
            daughtCalc[partTypes[0]+'_ET']+= daughtCalc[partTypes[2]+'_smeared']*missingFraction( daughtCalc['rdiff_'+partTypes[0]+partTypes[2]], daughtParams[partTypes[0]+'_L0Calo_HCAL_region'],clusters )
            daughtCalc[partTypes[2]+'_ET']+= daughtCalc[partTypes[0]+'_smeared']*missingFraction( daughtCalc['rdiff_'+partTypes[0]+partTypes[2]], daughtParams[partTypes[0]+'_L0Calo_HCAL_region'],clusters ) 
                                                                     


        Lc_ET_emulated_ = max(daughtCalc['pi_ET'],daughtCalc['p_ET'],daughtCalc['K_ET'])
        #print(i,Lc_ET_emulated_)
        if Lc_ET_emulated_ > 6100:
            Lc_ET_emulated_ = 6100
        if Lc_ET_emulated_ < 0:
            Lc_ET_emulated_ = 0

        Lc_ET_emulated_array[i] = Lc_ET_emulated_
        inputarray[i][0] = ev.Lc_PT
        inputarray[i][1] = ev.Lc_P
        inputarray[i][2] = daughtParams['K_L0Calo_HCAL_realET']
        inputarray[i][3] = daughtParams['pi_L0Calo_HCAL_realET']
        inputarray[i][4] = daughtParams['p_L0Calo_HCAL_realET']
        inputarray[i][5] = daughtCalc['rdiff_ppi']
        inputarray[i][6] = daughtCalc['rdiff_Kpi']
        inputarray[i][7] = daughtCalc['rdiff_Kp']

        if (i == 1):
          print(inputarray[1])


    bdt = joblib.load("bdt.joblib")
    y_predicted = bdt.predict(inputarray)

    Lc_L0Hadron_TOS_emulated = array('i',[0])
    br_Lc_L0Hadron_TOS = tout.Branch('Lc_L0Hadron_TOS_emulated', 
    Lc_L0Hadron_TOS_emulated, 'Lc_L0Hadron_TOS_emulated[1]/I')

    Lc_L0Hadron_ET_emulated = array('f',[0])
    br_Lc_L0Hadron_ET = tout.Branch('Lc_L0Hadron_ET_emulated', 
    Lc_L0Hadron_ET_emulated, 'Lc_L0Hadron_ET_emulated[1]/F')

    Lc_L0Hadron_BDT_correction = array('f',[0])
    br_Lc_L0Hadron_BDT_correction = tout.Branch('Lc_L0Hadron_BDT_correction', 
    Lc_L0Hadron_BDT_correction, 'Lc_L0Hadron_BDT_correction[1]/F')

    for i,ev in enumerate(tin):
        Lc_ET_BDT_corrected = Lc_ET_emulated_array[i]+y_predicted[i]
        Lc_L0Hadron_TOS_emulated[0] = 0
        Lc_L0Hadron_BDT_correction[0] = y_predicted[i]
        Lc_L0Hadron_ET_emulated[0] = Lc_ET_emulated_array[i]
        if Lc_ET_BDT_corrected > TOSthreshold:    
            Lc_L0Hadron_TOS_emulated[0] = 1
        tout.Fill()

    print("All events processed.")

    tout.Write()
    fout.Close()
    fin.Close()
