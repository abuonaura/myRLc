import sys, os
import ROOT as r
import pandas as pd
import numpy as np
from root_pandas import read_root
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo

filedir = '/disk/lhcb_data2/RLcMuonic2016/'

def ReadRootFile(dtype, polarity):
    df_list =[]
    varsON = ['Lb_L0Global_TIS','Lb_L0HadronDecision_TOS','Lc_Hlt1TrackMVADecision_TOS','Lc_Hlt1TwoTrackMVADecision_TOS','Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS','Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS','Lc_M','p_ProbNNp','p_ProbNNk','mu_PID*','*_P','*_PT','nTracks','runNumber','eventNumber','Lb_ISOLATION_*','mu_PX','mu_PY','mu_PZ','mu_ID','Lc_PX','Lc_PY','Lc_PZ']
    for df in read_root(filedir+'Data/Lb_'+dtype+'_'+polarity+'.root','tupleout/DecayTree', chunksize=100000, columns=varsON):
        df_list.append(df)
    return df_list

def L0TriggerData(df):
    check = lambda x,y: True if (x | y) else False
    df['L0'] = df.apply(lambda x: check(x['Lb_L0Global_TIS'],x['Lb_L0HadronDecision_TOS']),axis=1)
    return df

def HLT1TriggerData(df):
    check = lambda x,y: True if (x | y) else False
    df['HLT1'] = df.apply(lambda x: check(x['Lc_Hlt1TrackMVADecision_TOS'],x['Lc_Hlt1TwoTrackMVADecision_TOS']),axis=1)
    return df

def HLT2TriggerData(df,dtype):
    if dtype =='Data' or dtype == 'DataSS':
        df['HLT2'] = df.apply((lambda x: True if x['Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS'] else False),axis=1)
    if dtype =='FakeMu' or dtype=='FakeMuSS':
        df['HLT2'] = df.apply((lambda x: True if x['Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS'] else False),axis=1)
    return df

def TriggerData(df):
    check = lambda x,y,z: True if (x & y & z) else False
    df['Trigger'] = df.apply(lambda x: check(x['L0'],x['HLT1'],x['HLT2']),axis=1)
    return df

def LcMassCut(df):
    df['LcMass'] = df.apply(lambda x: True if 2230<= x['Lc_M'] <=2330 else False, axis=1)
    return df

def ApplyPIDCalibCuts(df):
    CheckNtrks = lambda x: True if x>0 and x<700 else False
    CheckPAndPT = lambda x,y: True if (x>0 and x<200000 and y>0 and y<60000) else False
    CheckAll = lambda x, y, z, v, t: True if x & y & z & v & t else False
    particles = ['mu','p','pi','K']
    df['PIDCalib'] = df.apply(lambda x:True if (CheckNtrks(x['nTracks'])& CheckPAndPT(x['mu_P'],x['mu_PT']) &
                              CheckPAndPT(x['pi_P'],x['pi_PT']) & CheckPAndPT(x['p_P'],x['p_PT']) &
                              CheckPAndPT(x['K_P'],x['K_PT']) )else False, axis=1)
    return df

def ApplyMuCuts(df,dtype):
    if dtype=='Data'or dtype=='DataSS':
        df['MuCuts'] = df.apply(lambda x: True if (x['mu_PIDmu']>2 and x['mu_PIDmu']-x['mu_PIDK']>2 and
                                                  x['mu_PIDmu']-x['mu_PIDp']>2) else False, axis=1)
    else:
        df['MuCuts'] = True
    return df

def GetFinalPreselection(df):
    df['Preselection'] = df.apply(lambda x: True if x['Trigger']& x['LcMass'] &x['PIDCalib'] & 
                                  x['MuCuts'] else False,axis=1)
    return df


def LoadBDTdf(dtype,polarity):
    df_list_bdt=[]
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
    for df_bdt in read_root(bdtfname,'DecayTree', chunksize=100000):
        df_list_bdt.append(df_bdt)
    return df_list_bdt

def MergeDataFrames(df,df_bdt):
    mergedDf = df.merge(df_bdt['bdt'],left_index=True, right_index=True)
    return mergedDf

BDTcut=0.7
def PassBDT(df,BDTcut):
    df['PassBDT'] = df.apply(lambda x: True if x['bdt']>BDTcut else False, axis=1)
    return df

def RemoveDDstar(df):
    df['NoDDstar'] = df.apply(lambda x: True if x['p_ProbNNp']- x['p_ProbNNk']>0 else False,axis=1)
    return df

def FinalSelection(df):
    df['FinalSel']= df.apply(lambda x: True if x['Preselection']&x['PassBDT']&x['NoDDstar'] else False, axis=1)
    return df

ISOBDTcut =0.35
ISOBDT2cut=0.2

#Masses pi, K, p, mu, Lc
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)

def CheckIfIsKenriched(df,n):
    for a in range (len(df)):
        i = a+n*100000
        #print(i)
        BDT = df.loc[i,'Lb_ISOLATION_BDT']
        BDT2 = df.loc[i,'Lb_ISOLATION_BDT2']
        PIDp = df.loc[i,'Lb_ISOLATION_PIDp']
        PIDK = df.loc[i,'Lb_ISOLATION_PIDK']
        Ch = df.loc[i,'Lb_ISOLATION_CHARGE']
        PIDp2 = df.loc[i,'Lb_ISOLATION_PIDp2']
        PIDK2 =df.loc[i,'Lb_ISOLATION_PIDK2']
        Ch2 = df.loc[i,'Lb_ISOLATION_CHARGE2']
        muCh = -df.loc[i,'mu_ID']/13
        isK1, isK2=0,0

        if BDT>ISOBDTcut and BDT2>ISOBDT2cut:
            m1 = m_pi
            m2 = m_pi
            #I retrieve the momenta of the 2 anti-isolated particles
            #------------------------
            #I will do in the following several mass hypothesis for these 2 particles to make
            #the Kenriched sample as cleaner as possible
            #------------------------
            #I want to discard events coming from LcStar ->Lcpipi
            #-> I assume that the 2 anti-isolated particles are 2 pions
            p1=np.array(df.loc[i,['Lb_ISOLATION_PX','Lb_ISOLATION_PY','Lb_ISOLATION_PZ']])
            E1pi =np.sqrt(p1.dot(p1.transpose())+m_pi*m_pi)
            p2=np.array(df.loc[i,['Lb_ISOLATION_PX2','Lb_ISOLATION_PY2','Lb_ISOLATION_PZ2']])
            E2pi =np.sqrt(p2.dot(p2.transpose())+m_pi*m_pi)
            #print(E1pi,E2pi)
            pLc=np.array(df.loc[i,['Lc_PX','Lc_PY','Lc_PZ']])
            ELc = np.sqrt(pLc.dot(pLc.transpose())+m_Lc*m_Lc)
            ELc12 = E1pi+E2pi+ELc
            pLc12 = p1+p2+pLc
            mLc12= np.sqrt(ELc12*ELc12 - pLc12.dot(pLc12.transpose()))
            if PIDp>4.:
                if (Ch==-muCh or (Ch==muCh and (PIDp-PIDK)<0.)):
                    m1 = m_K
                    isK1=1
                if (Ch==muCh and PIDp-PIDK>0.):#I assume it is a pion
                    m1 = m_pi
            if PIDK2>4.:
                if (Ch2==-muCh or (Ch2==muCh and (PIDp2 - PIDK2)<0.)):
                    #I also assume this particle is a K
                    m2 = m_K
                    isK2=1
                if (Ch2==muCh and PIDp2-PIDK2>0.):
                    m2 = m_pi
            #-> I compute the energy of the 2 particles:
            E1 = np.sqrt(p1.dot(p1.transpose()) + m1*m1)
            E2 = np.sqrt(p2.dot(p2.transpose()) + m2*m2)
            #-> I retrieve the momentum of the muon
            pmu = np.array(df.loc[i,['mu_PX','mu_PY','mu_PZ']])
            Emu = np.sqrt(pmu.dot(pmu.transpose()) + m_mu*m_mu)
            #-> I compute the total momentum and energy of the 4 particles:
            pTOT = p1 + p2 + pmu + pLc
            ETOT = E1 + E2 + Emu + ELc
            #-> I evaluate the invariant mass
            mTOT = np.sqrt(ETOT*ETOT - pTOT.dot(pTOT.transpose()))

            if mLc12>2770 and mTOT<5620 and (isK1==1 or isK2==1):
                df.at[i,'isKenriched'] = True
                #print(mLc12[i])
            else:
                df.at[i,'isKenriched'] = False
        else:
            df.at[i,'isKenriched'] = False
    return df

def CheckIfIsIsolated(df):
    df['isIsolated'] = df.apply(lambda x: True if x['Lb_ISOLATION_BDT']<ISOBDTcut
                                else False, axis=1)
    return df



dflist = ReadRootFile('DataSS','MagUp')
dflist_bdt = LoadBDTdf('DataSS','MagUp')
dflist_final = []
for n, df in enumerate(dflist):
    print(n)
    if os.path.isfile('dataSS_'+str(n)+'.h5'):
        df = pd.read_hdf('dataSS_'+str(n)+'.h5', 'df')
    else:
        df.to_hdf('dataSS_'+str(n)+'.h5',key='df',mode='w')
    df_bdt = dflist_bdt[n]
    print('Applying L0 Trigger: ')
    df = L0TriggerData(df)
    print('Applying HLT1 Trigger: ')
    df = HLT1TriggerData(df)
    print('Applying HLT2 Trigger: ')
    df = HLT2TriggerData(df,'DataSS')
    print('Applying Trigger: ')
    df = TriggerData(df)
    print('Apply Lc Mass cut: ')
    df = LcMassCut(df)
    print('Apply PIDcalib cut: ')
    df = ApplyPIDCalibCuts(df)
    print('Apply muPID cut: ')
    df= ApplyMuCuts(df,'DataSS')
    print('Apply final preselection: ')
    df=GetFinalPreselection(df)
    print('Adding BDT: ')
    mergedDf = df.merge(df_bdt['bdt'],left_index=True, right_index=True)
    print('Cut on BDT: ')
    mergedDf=PassBDT(mergedDf,BDTcut)
    print('DDstar cut: ')
    mergedDf= RemoveDDstar(mergedDf)
    print('Final selection')
    mergedDf = FinalSelection(mergedDf)
    print('IsKenriched')
    mergedDf = CheckIfIsKenriched(mergedDf,n)
    print('IsISolated')
    mergedDf = CheckIfIsIsolated(mergedDf)
    mergedDf.head()
    dflist_final.append(mergedDf)

print('concatenate')
finalDf = pd.concat(dflist_final)
finalDf.to_root('Lb_DataSS_Df_preselection.root','DecayTree')


