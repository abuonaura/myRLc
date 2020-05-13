import ROOT as r
import math as m
import sys, os

from sklearn import datasets
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.metrics import classification_report, roc_auc_score
#from sklearn.externals import joblib
import joblib
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
from array import array
import time
import pandas as pd


function_code = '''

int mombins = 10;
int ptbins = 6;
double pt_low = 0;
double pt_high = 15000;
double mom_low = 0;
double mom_high = 1e5;
TFile* momhistsfile = new TFile("L0Hadron_TOS/momhists.root");
std::vector< std::vector< TH1D* > > momPThists;

TRandom3 rando;

TFile* fcluster_hists = new TFile("L0Hadron_TOS/two_particle_HCAL_cluster_hists.root");
TH1D* hshared_inner = (TH1D*)fcluster_hists->Get("shared_with_radial_inner");
TH1D* hshared_outer = (TH1D*)fcluster_hists->Get("shared_with_radial_outer");
TH1D* hmissing_inner = (TH1D*)fcluster_hists->Get("missing_with_radial_inner");
TH1D* hmissing_outer = (TH1D*)fcluster_hists->Get("missing_with_radial_outer");

int init(){
  for (int i=0; i<mombins; i++){
    std::vector < TH1D* > PThists; 
    for (int j=0; j<mombins; j++){
      TString hname = "hdiff_";
      hname+=i;
      hname+="_";
      hname+=j;
      PThists.push_back((TH1D*)momhistsfile->Get(hname));
    }
    momPThists.push_back(PThists);
  }
 return true;
}

double random_smearing(double P,double PT,double realET){
  int Pbin = int(floor(P/((mom_high-mom_low)/mombins)));       
  int PTbin = int(floor(PT/((pt_high-pt_low)/ptbins)));
  if (Pbin < 0) Pbin=0;
  if (PTbin < 0) PTbin=0;
  if (Pbin >= mombins) Pbin=mombins-1;
  if (PTbin >= ptbins) PTbin=ptbins-1;
  TH1D* momPThist = momPThists[Pbin][PTbin];
  double rando = momPThist->GetRandom();
  double smearedET = realET*(1-rando);
  if (smearedET < 0) smearedET=0;
  if (smearedET > 6100) smearedET=6100;
  return smearedET;
  }

double rdiff(double xproj1,double yproj1,double xproj2,double yproj2){
  return sqrt((xproj1-xproj2)*(xproj1-xproj2)+(yproj1-yproj2)*(yproj1-yproj2));
}

int isShared(double rdiff, int region1, int region2){
  if (region1!=region2) return false;
  TH1D* hshared = hshared_outer;
  if (region1 == 1) hshared = hshared_inner;
  double fracshared = hshared->GetBinContent(hshared->FindBin(rdiff));
  if (rando.Uniform() < fracshared){
   return 1;
  }
  else return 0;
}

double missingFraction(double rdiff,int region1,int region2){
  if (region1!=region2) return false;
  TH1D* hmissing = hmissing_outer;
  if (region1 == 1) hmissing = hmissing_inner;
  double fracmissing = hmissing->GetBinContent(hmissing->FindBin(rdiff))/hmissing->GetBinContent(10);
  return fracmissing;
}

double calcET(double smear1,double smear2,double smear3,int shared12,int shared13,double missing12,double missing13){
  double ET = smear1;
  if (shared12) ET = (smear1+smear2)*missing12;
  if (shared13) ET+= smear3*missing13;
  return ET;
}

double calcLcET(double KET,double pET, double piET){
  if ((KET >= piET) && (KET >= pET)) return KET;
  if ((pET >= KET) && (pET >= piET)) return pET;
  if ((piET >= pET) && (piET >= KET)) return piET;
}

int calcTOS(double KET, double pET, double piET){ 
  double TOSthreshold = 3744.0;
  if ((KET > TOSthreshold) || (pET > TOSthreshold) || (piET > TOSthreshold)){
    return 1;
  }
  else {
    return 0;
  }
}
'''




def main(MC_fileName = "../../../fromSimone/L0HLT1emulationVars/Lb_Lcmunu_MagUp.root"):

 





    rando = r.TRandom3()
    fcluster_hists = r.TFile.Open("L0Hadron_TOS/two_particle_HCAL_cluster_hists.root")
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
    tin.SetBranchStatus("Lc_P",1)
    tin.SetBranchStatus("Lc_PT",1)
    r.gInterpreter.Declare(function_code)
    print(r.init())
    df0 = r.RDataFrame(tin)


    start = time.time()
    df1 = df0.Define("K_smeared","random_smearing(K_P,K_PT,K_L0Calo_HCAL_realET)") 
    df2 = df1.Define("pi_smeared","random_smearing(pi_P,pi_PT,pi_L0Calo_HCAL_realET)") 
    df3 = df2.Define("p_smeared","random_smearing(p_P,p_PT,p_L0Calo_HCAL_realET)") 
    df4 = df3.Define("rdiffKpi","rdiff(K_L0Calo_HCAL_xProjection,K_L0Calo_HCAL_yProjection,pi_L0Calo_HCAL_xProjection,pi_L0Calo_HCAL_yProjection)") 
    df5 = df4.Define("rdiffKp","rdiff(K_L0Calo_HCAL_xProjection,K_L0Calo_HCAL_yProjection,p_L0Calo_HCAL_xProjection,p_L0Calo_HCAL_yProjection)") 
    df6 = df5.Define("rdiffpip","rdiff(pi_L0Calo_HCAL_xProjection,pi_L0Calo_HCAL_yProjection,p_L0Calo_HCAL_xProjection,p_L0Calo_HCAL_yProjection)")
    #df7 = df6.Define("sharedKpi","random_smearing(rdiffKpi,K_L0Calo_HCAL_region,pi_L0Calo_HCAL_region)")
    df7 = df6.Define("sharedKpi","isShared(rdiffKpi,K_L0Calo_HCAL_region,pi_L0Calo_HCAL_region)")
    df8 = df7.Define("sharedpip","isShared(rdiffpip,pi_L0Calo_HCAL_region,p_L0Calo_HCAL_region)")
    df9 = df8.Define("sharedKp","isShared(rdiffKp,K_L0Calo_HCAL_region,p_L0Calo_HCAL_region)")
    df10 = df9.Define("missingKpi","missingFraction(rdiffKpi,K_L0Calo_HCAL_region,pi_L0Calo_HCAL_region)")
    df11 = df10.Define("missingpip","missingFraction(rdiffpip,pi_L0Calo_HCAL_region,p_L0Calo_HCAL_region)")
    df12 = df11.Define("missingKp","missingFraction(rdiffKp,K_L0Calo_HCAL_region,p_L0Calo_HCAL_region)")
    df13 = df12.Define("KET","calcET(K_smeared,p_smeared,pi_smeared,sharedKp,sharedKpi,missingKp,missingKpi)")
    df14 = df13.Define("pET","calcET(p_smeared,K_smeared,pi_smeared,sharedKp,sharedpip,missingKp,missingpip)")
    df15 = df14.Define("piET","calcET(pi_smeared,p_smeared,K_smeared,sharedpip,sharedKpi,missingpip,missingKpi)")
    df16 = df15.Define("LcET","calcLcET(piET,pET,KET)")
    #df7 = df6.Define("sharedKpi","rdiff(pi_L0Calo_HCAL_xProjection,pi_L0Calo_HCAL_yProjection,p_L0Calo_HCAL_xProjection,p_L0Calo_HCAL_yProjection)")
    print("Writing the new tree ...")

    column_name_vector = r.std.vector('string')()
    #column_name_vector.push_back("K_smeared")
    #column_name_vector.push_back("pi_smeared")
    #column_name_vector.push_back("p_smeared")
    column_name_vector.push_back("Lc_PT")
    column_name_vector.push_back("Lc_P")
    column_name_vector.push_back("K_L0Calo_HCAL_realET")
    column_name_vector.push_back("pi_L0Calo_HCAL_realET")
    column_name_vector.push_back("p_L0Calo_HCAL_realET")
    column_name_vector.push_back("rdiffKpi")
    column_name_vector.push_back("rdiffKp")
    column_name_vector.push_back("rdiffpip")
    column_name_vector.push_back("LcET")

    #df16.Snapshot("DecayTree",'_wHLT1Emulation.root', column_name_vector)
    #print("Tree written.")
    bdt_vars = ["Lc_PT","Lc_P","K_L0Calo_HCAL_realET","pi_L0Calo_HCAL_realET","p_L0Calo_HCAL_realET","rdiffpip","rdiffKpi","rdiffKp"]
    bdtinputdict = df16.AsNumpy(columns=bdt_vars)
    LcETraw = df16.AsNumpy(columns=["LcET"])["LcET"]
 
    inputarray = bdtinputdict["Lc_PT"]
    bdt_vars.remove("Lc_PT")
    for var in bdt_vars:
      inputarray = np.c_[inputarray,bdtinputdict[var]]


    
  

    '''Lc_ET_emulated_array[i] = Lc_ET_emulated_
    inputarray[i][0] = npy2["Lc_PT"]
    inputarray[i][1] = ev.Lc_P
    inputarray[i][2] = daughtParams['K_L0Calo_HCAL_realET']
    inputarray[i][3] = daughtParams['pi_L0Calo_HCAL_realET']
    inputarray[i][4] = daughtParams['p_L0Calo_HCAL_realET']
    inputarray[i][5] = daughtCalc['rdiff_ppi']
    inputarray[i][6] = daughtCalc['rdiff_Kpi']
    inputarray[i][7] = daughtCalc['rdiff_Kp']'''

    bdt = joblib.load("L0Hadron_TOS/bdt.joblib")
    y_predicted = bdt.predict(inputarray)

    corrected_TOS = np.add(LcETraw,y_predicted)
 
    TOSthreshold = 3744.0
    corrected_TOS[corrected_TOS < TOSthreshold] = 0
    corrected_TOS[corrected_TOS > TOSthreshold] = 1
  
    corrected_TOS = corrected_TOS.astype(int)
    

    corrected_TOS.dtype = [('Lc_L0Hadron_TOS_emulated', 'int')]

    from root_numpy import array2root



    array2root(corrected_TOS,MC_fileName[0:-5]+"_wL0TOSEmulation.root", "DecayTree", mode='recreate')

 
    print('time taken for whole tree',time.time()-start)


