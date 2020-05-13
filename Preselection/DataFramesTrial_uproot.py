#!/usr/bin/env python3
import uproot
import sys, os
import pandas as pd
import numpy as np
from root_pandas import to_root
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo
from AddSweights import *

filedir = '/disk/lhcb_data2/RLcMuonic2016/'


def GetDataframes(dtype,polarity):
    varsON = ['Lb_L0Global_TIS','Lb_L0HadronDecision_TOS','Lc_Hlt1TrackMVADecision_TOS',
              'Lc_Hlt1TwoTrackMVADecision_TOS','Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS',
              'Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS','Lc_M','p_ProbNNp','p_ProbNNk','mu_PID*','*_P',
              '*_PT','nTracks','runNumber','eventNumber','Lb_ISOLATION_*','mu_PX','mu_PY','mu_PZ','mu_ID',
              'Lc_PX','Lc_PY','Lc_PZ']
    df_list=[]
    for df in uproot.pandas.iterate(filedir+'Data/Lb_'+dtype+'_'+polarity+'.root','tupleout/DecayTree',
                                varsON, entrysteps=500000):
        df_list.append(df)
    return df_list

def L0TriggerData(df):
    df['L0'] = False
    df.loc[df.Lb_L0Global_TIS | df.Lb_L0HadronDecision_TOS,'L0'] = True
    return df

def HLT1TriggerData(df):
    df['HLT1'] = False
    df.loc[df.Lc_Hlt1TrackMVADecision_TOS | df.Lc_Hlt1TwoTrackMVADecision_TOS,'HLT1'] = True
    return df

def HLT2TriggerData(df,dtype):
    df['HLT2'] = False
    if dtype =='Data' or dtype == 'DataSS':
        df.loc[df.Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS,'HLT2'] = True
    if dtype =='FakeMu' or dtype=='FakeMuSS':
        df.loc[df.Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS,'HLT2'] = True
    return df

def TriggerData(df):
    df['Trigger'] = False
    df.loc[df.L0 & df.HLT1 & df.HLT2,'Trigger'] = True
    return df

def LcMassCut(df):
    df['LcMass'] = False
    df.loc[(df.Lc_M>2230) & (df.Lc_M<2330),'LcMass'] = True
    return df

def ApplyPIDCalibCuts(df):
    df['PIDCalib'] = False
    df.loc[df.nTracks.between(0,700)&df.mu_P.between(0,200000)&df.mu_PT.between(0,60000)&
          df.pi_P.between(0,200000)&df.pi_PT.between(0,60000)&
          df.p_P.between(0,200000)&df.p_PT.between(0,60000)&
          df.K_P.between(0,200000)&df.K_PT.between(0,60000),'PIDCalib'] = True
    return df

def ApplyMuCuts(df,dtype):
    if dtype=='Data'or dtype=='DataSS':
        df['MuCuts'] = False
        df.loc[(df.mu_PIDmu>2) & ((df.mu_PIDmu -df.mu_PIDK)>2) & ((df.mu_PIDmu -df.mu_PIDp)>2) & (df.mu_PIDe<1),'MuCuts']=True
    else:
        df['MuCuts'] = True
    return df

def GetFinalPreselection(df):
    df['Preselection'] = False
    df.loc[df.Trigger & df.LcMass & df.MuCuts & df.PIDCalib,'Preselection'] = True
    return df

def LoadBDTdf(dtype,polarity):
    dflist_bdt=[]
    ifname = filedir+'Data/Lb_'+dtype+'_'+polarity+'.root'
    bdtfname = ifname[0:-5]+'_MVA.root'
    if os.path.isfile(bdtfname):
        print('BDT file already created')
    else:
        print()
        print('>>>   Creating file with BDT variable')
        print()
        AddBDTinfo(ifname, 'tupleout/DecayTree', bdtfname, 'Data',
                   pickled_model_path = '../PIDGen_PIDCalib_MVA/xgb_reg.pkl')
    for df in uproot.pandas.iterate(bdtfname,'DecayTree',
                            'bdt', entrysteps=500000):
        dflist_bdt.append(df)
    return dflist_bdt

BDTcut=0.7
def PassBDT(df,BDTcut):
    df['PassBDT'] = False
    df.loc[df.bdt>BDTcut,'PassBDT']=True
    return df

def RemoveDDstar(df):
    df['NoDDstar'] = False
    df.loc[((df.p_ProbNNp - df.p_ProbNNk)>0),'NoDDstar']=True
    return df

def FinalSelection(df):
    df['FinalSel']= False
    df.loc[df.Preselection & df.PassBDT & df.NoDDstar,'FinalSel']=True
    return df

ISOBDTcut =0.35
ISOBDT2cut=0.2

#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)

def CheckIfIsKenriched(df):
    df1 = df[(df.Lb_ISOLATION_BDT>ISOBDTcut) & (df.Lb_ISOLATION_BDT2>ISOBDT2cut)].copy()
    df1['E1pi'] = np.sqrt(df1["Lb_ISOLATION_PX"]**2 + df1['Lb_ISOLATION_PY']**2 + df1['Lb_ISOLATION_PZ']**2 + m_pi**2)
    df1['E2pi'] = np.sqrt(df1["Lb_ISOLATION_PX2"]**2 + df1['Lb_ISOLATION_PY2']**2 + df1['Lb_ISOLATION_PZ2']**2 + m_pi**2)
    df1['ELc'] =np.sqrt(df1['Lc_PX']**2+df1['Lc_PY']**2+df1['Lc_PZ']**2 + m_Lc**2)
    df1['pLc12_x'] = df1['Lc_PX']+df1["Lb_ISOLATION_PX"]+df1["Lb_ISOLATION_PX2"]
    df1['pLc12_y'] = df1['Lc_PY']+df1["Lb_ISOLATION_PY"]+df1["Lb_ISOLATION_PY2"]
    df1['pLc12_z'] = df1['Lc_PZ']+df1["Lb_ISOLATION_PZ"]+df1["Lb_ISOLATION_PZ2"]
    df1['mLc12'] = np.sqrt((df1['ELc']+df1['E1pi']+df1['E2pi'])**2 -
                           (df1['pLc12_x']**2+df1['pLc12_y']**2+df1['pLc12_z']**2))
    df1['muCharge'] = -df1['mu_ID']/13
    df1['PIDdiff'] =df1['Lb_ISOLATION_PIDp'] - df1['Lb_ISOLATION_PIDK']
    df1['PIDdiff2'] =df1['Lb_ISOLATION_PIDp2'] - df1['Lb_ISOLATION_PIDK2']
    df1['m1']= m_pi
    df1['m2']= m_pi
    df1.loc[(df1.Lb_ISOLATION_PIDK>4.)&((df1.Lb_ISOLATION_CHARGE==-df1.muCharge)|
                                        ((df1.Lb_ISOLATION_CHARGE==df1.muCharge) 
                                         & (df1.PIDdiff<0))),'m1']= m_K
    df1.loc[(df1.Lb_ISOLATION_PIDK2>4.)&((df1.Lb_ISOLATION_CHARGE2==-df1.muCharge)|
                                        ((df1.Lb_ISOLATION_CHARGE2==df1.muCharge) 
                                         & (df1.PIDdiff2<0))),'m2']= m_K
    df1['E1'] = np.sqrt(df1["Lb_ISOLATION_PX"]**2 + df1['Lb_ISOLATION_PY']**2 + df1['Lb_ISOLATION_PZ']**2 +df1['m1']**2)
    df1['E2'] = np.sqrt(df1["Lb_ISOLATION_PX2"]**2 + df1['Lb_ISOLATION_PY2']**2 + df1['Lb_ISOLATION_PZ2']**2 +df1['m2']**2)
    df1['Emu'] = np.sqrt(df1.mu_PX**2 + df1.mu_PY**2 + df1.mu_PZ**2 +m_mu**2)
    df1['ETOT'] = df1['ELc']+df1['E1']+df1['E2']+df1['Emu']
    df1['pTOT_x'] = df1['pLc12_x']+df1.mu_PX
    df1['pTOT_y'] = df1['pLc12_y']+df1.mu_PY
    df1['pTOT_z'] = df1['pLc12_z']+df1.mu_PZ
    df1['mTOT'] = np.sqrt(df1.ETOT**2 -(df1.pTOT_x**2+df1.pTOT_y**2+df1.pTOT_z**2))
    df1['isKenriched']=False
    #print(df1['isKenriched'].dtype)
    df1.loc[(df1.mLc12>2700) & (df1.mTOT<5620) & ((df1.m1==m_K) | (df1.m2==m_K)), 'isKenriched']=True 
    #print(df1['isKenriched'].dtype)
    df1=df1.drop(df1.loc[:, 'E1pi':'mTOT'].columns,axis=1)
    #print(df1['isKenriched'].dtype)
    dfT = pd.merge(df,df1['isKenriched'],how='left',left_index=True, right_index=True)
    #print(dfT['isKenriched'].dtype)
    dfT.loc[dfT['isKenriched'].isnull(),'isKenriched']=False
    dfT['isKenriched'] = dfT['isKenriched'].astype('bool')
    return dfT

def CheckIfIsIsolated(df):
    df['isIsolated'] = False
    df.loc[df.Lb_ISOLATION_BDT<ISOBDTcut,'isIsolated']=True
    return df


def CreatePreselectionTree(dtype,polarity):
    df_list = GetDataframes(dtype,polarity)
    #df = pd.read_csv("df0.csv")
    dflist_final = []
    dflist_final_Full = []
    dflist_final_Iso = []
    dflist_final_Kenr = []
    dflist_bdt = LoadBDTdf(dtype,polarity)
    print('Total number of dataframes to analyse: ', len(df_list))
    for n, df in enumerate(df_list):
    #for n in range(0,1):
        if n%10==0:
            print(n)
        df_bdt = dflist_bdt[n]
        #print('Applying L0 Trigger: ')
        df = L0TriggerData(df)
        #print('Applying HLT1 Trigger: ')
        df = HLT1TriggerData(df)
        #print('Applying HLT2 Trigger: ')
        df = HLT2TriggerData(df,dtype)
        #print('Applying Trigger: ')
        df = TriggerData(df)
        #print('Apply Lc Mass cut: ')
        df = LcMassCut(df)
        #print('Apply PIDcalib cut: ')
        df = ApplyPIDCalibCuts(df)
        #print('Apply muPID cut: ')
        df= ApplyMuCuts(df,dtype)
        #print('Apply final preselection: ')
        df=GetFinalPreselection(df)
        #print('Adding BDT: ')
        df2 = df.merge(df_bdt['bdt'],left_index=True, right_index=True)
        #print('Cut on BDT: ')
        df2 = PassBDT(df2,BDTcut)
        #print('DDstar cut: ')
        df2 = RemoveDDstar(df2)
        #print('Final selection')
        df2 = FinalSelection(df2)
        #print('IsKenriched')
        dfT = CheckIfIsKenriched(df2)
        #print(dfT['isKenriched'].dtype)
        #print('IsISolated')
        dfT = CheckIfIsIsolated(dfT)

       
        #Create Kenr dataframe
        dfKenr =  dfT.copy()
        dfKenr = dfKenr.loc[(dfKenr.FinalSel==True) & (dfKenr.isIsolated==False) & (dfKenr.isKenriched==True)]
        #dfKenr = dfKenr.drop(dfKenr.loc[:,'L0':'isIsolated'].columns,axis=1)
        dflist_final_Kenr.append(dfKenr)
        
        #Create Iso dataframe
        dfIso =  dfT.copy()
        dfIso = dfIso.loc[(dfIso.FinalSel==True) & (dfIso.isIsolated==True) & (dfIso.isKenriched==False)]
        #dfIso = dfIso.drop(dfIso.loc[:,'L0':'isIsolated'].columns,axis=1)
        dflist_final_Iso.append(dfIso)
        
        #create Full dataframe
        dfFull =  dfT.copy()
        dfFull = dfFull.loc[dfFull.FinalSel==True]
        #dfFull = dfFull.drop(dfFull.loc[:,'L0':'isIsolated'].columns,axis=1)
        dflist_final_Full.append(dfFull)
        
        #dfT = dfT.drop(dfT.loc[:,'Lb_L0Global_TIS':'Lc_PZ'].columns,axis=1)
        dflist_final.append(dfT)

    #print('concatenate')
    finalDf = pd.concat(dflist_final)
    #finalDf.to_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Df_preselection.root','DecayTree')
    #Create Full root file
    finalDf_Full = pd.concat(dflist_final_Full)
    finalDf_Full.to_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_full.root','DecayTree')
    ofname = ComputeSweights(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_full.root',dtype,polarity,'full')
    print('Created: '+ofname)
    f_sw = uproot.open(ofname)
    t_sw = f_sw['DecayTree']
    df_Full_sw = t_sw.pandas.df()
    df_Full_sw = df_Full_sw.set_index('__index__')
    #Create Iso root file
    finalDf_Iso = pd.concat(dflist_final_Iso)
    finalDf_Iso.to_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_iso.root','DecayTree')
    ofname = ComputeSweights(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_iso.root',dtype,polarity,'iso')
    print('Created: '+ofname)
    f_sw = uproot.open(ofname)
    t_sw = f_sw['DecayTree']
    df_Iso_sw = t_sw.pandas.df()
    df_Iso_sw = df_Iso_sw.set_index('__index__')
    #Create Kenr root file
    finalDf_Kenr = pd.concat(dflist_final_Kenr)
    finalDf_Kenr.to_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_Kenr.root','DecayTree')
    ofname = ComputeSweights(filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_Kenr.root',dtype,polarity,'Kenriched')
    print('Created: '+ofname)
    f_sw = uproot.open(ofname)
    t_sw = f_sw['DecayTree']
    df_Kenr_sw = t_sw.pandas.df()
    df_Kenr_sw = df_Kenr_sw.set_index('__index__')
    #merge all the results together
    finalDf_sw = finalDf.merge(df_Full_sw[['sw_sig','sw_bkg']],how='left',left_index=True,right_index=True)
    finalDf_sw = finalDf_sw.rename(columns={'sw_sig':'sw_sig_Full','sw_bkg':'sw_bkg_Full'})
    finalDf_sw = finalDf_sw.merge(df_Iso_sw[['sw_sig','sw_bkg']],how='left',left_index=True,right_index=True)
    finalDf_sw = finalDf_sw.rename(columns={'sw_sig':'sw_sig_Iso','sw_bkg':'sw_bkg_Iso'})
    finalDf_sw = finalDf_sw.merge(df_Kenr_sw[['sw_sig','sw_bkg']],how='left',left_index=True,right_index=True)
    finalDf_sw = finalDf_sw.rename(columns={'sw_sig':'sw_sig_Kenr','sw_bkg':'sw_bkg_Kenr'})
    finalDf_sw = finalDf_sw.drop(finalDf_sw.loc[:,'Lb_L0Global_TIS':'Lc_PZ'].columns,axis=1)

    finalDf_sw.to_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Df_preselection_sw.root','DecayTree')

    for sm in ['full','iso','Kenr']:
        os.system('rm '+filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_'+sm+'.root')
        os.system('rm '+filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected_'+sm+'_sw.root')


    return finalDf,finalDf_sw, finalDf_Full, finalDf_Iso, finalDf_Kenr

def CreateFinalTree(dtype,polarity,sample):
    f = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'.root','READ')
    t = f.Get('tupleout/DecayTree')
    fsel = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Df_preselection.root','READ')
    tsel = fsel.Get('DecayTree')
    t.AddFriend(tsel)
    ofname = filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected'+suffix[sample]
    of = r.TFile(ofname,'RECREATE')
    ot = t.CloneTree()
    sw_bkg = np.zeros(1, dtype=float)
    sw_sig = np.zeros(1, dtype=float)
    sw_bkg_b = ot.Branch("sw_bkg",sw_bkg,"sw_bkg/D")
    sw_sig_b = ot.Branch("sw_sig",sw_sig,"sw_sig/D")
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if sample=='full':
            if tsel.FinalSel==1:
                sw_sig[0] = tsel.sw_sig_Full
                sw_bkg[0] = tsel.sw_bkg_Full
                ot.Fill()
        if sample=='iso':
            if tsel.FinalSel==1 and tsel.isIsolated==1 and tsel.isKenriched==0:
                sw_sig[0] = tsel.sw_sig_Iso
                sw_bkg[0] = tsel.sw_bkg_Iso
                ot.Fill()
        if sample=='Kenriched':
            if tsel.FinalSel==1 and tsel.isIsolated==0 and tsel.isKenriched==1:
                sw_sig[0] = tsel.sw_sig_Kenr
                sw_bkg[0] = tsel.sw_bkg_Kenr
                ot.Fill()

    of.Write()
    of.Close()
    print('Created FINAL TREE: '+filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected'+suffix[sample])
    return


if __name__ == "__main__":
    dtype = ['Data','DataSS','FakeMu','FakeMuSS']
    #dtype = ['DataSS']
    polarities=['MagUp','MagDown']
    #polarities=['MagDown']
    samples = ['full','iso','Kenriched']
    suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}
    for polarity in polarities:
        print(polarity)
        for dt in dtype:
            print('   ', dt)
            #CreatePreselectionTree(dt,polarity)
            for sample in samples:
                print('      ',sample)
                CreateFinalTree(dt,polarity,sample)

