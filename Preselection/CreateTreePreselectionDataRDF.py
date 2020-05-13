import sys, os
import pandas as pd
import numpy as np
from root_pandas import to_root
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo
from AddSweights import *

filedir = '/disk/lhcb_data2/RLcMuonic2016/DataRDF/'

def GetBDTfile(dtype,polarity):
    ifname = filedir+'Lb_'+dtype+'_'+polarity+'.root'
    bdtfname = ifname[0:-5]+'_MVA.root'
    if os.path.isfile(bdtfname):
        print('BDT file already created')
    else:
        print()
        print('>>>   Creating file with BDT variable')
        print()
        AddBDTinfo(ifname, 'tupleout/DecayTree', bdtfname, 'Data',
                   pickled_model_path = '../PIDGen_PIDCalib_MVA/xgb_reg.pkl')
    return bdtfname

func_code = '''
bool L0Trigger(bool L0TOS, bool L0TIS)
{
    if(L0TOS || L0TIS)
        return true;
    else
        return false;
}

bool HLT1Trigger(bool HLT1_1, bool HLT1_2)
{
    if(HLT1_1 || HLT1_2)
        return true;
    else
        return false;
}

bool HLT2Trigger(double HLT2_TOS)
{
    if(HLT2_TOS)
        return true;
    else
        return false;
}

bool PassTrigger(bool L0, bool HLT1, bool HLT2)
{
    if(L0 && HLT1 && HLT2)
        return true;
    else
        return false;
}

bool LcMassCut(double m)
{
    if(m>2230 && m<2330)
        return true;
    else
        return false;
}

bool ApplyPIDCalibCuts(double nTracks, double K_P, double K_PT, double p_P, double p_PT,double pi_P, double pi_PT, double mu_P, double mu_PT)
{
    if(nTracks>0 && nTracks<700)
    {
        if(K_P>0 && K_P<2.E5 && K_PT>0 && K_PT<6E4)
        {
            if(p_P>0 && p_P<2.E5 && p_PT>0 && p_PT<6E4)
            {
                if(pi_P>0 && pi_P<2.E5 && pi_PT>0 && pi_PT<6E4)
                {
                    if(mu_P>0 && mu_P<2.E5 && mu_PT>0 && mu_PT<6E4)
                        return true;
                }
            }
        }
    }
    return false;
}

bool ApplyMuCuts(double mu_PIDmu, double mu_PIDK, double mu_PIDp, double mu_PIDe)
{
    if(mu_PIDmu>2. && (mu_PIDmu- mu_PIDK)>2 && (mu_PIDmu - mu_PIDp)>2. && mu_PIDe<1)
        return true;
    else
        return false;
}

bool GetFinalPreselection(bool Trigger, bool PIDCalibCuts, bool LcMass, bool MuCuts)
{
    if (Trigger && PIDCalibCuts && LcMass && MuCuts)
        return true;
    else
        return false;
}

bool PassBDT(double bdt)
{
    double BDTcut = 0.7;
    if(bdt>BDTcut)
        return true;
    else
        return false;
}

bool RemoveDDstar(double p_ProbNNp, double p_ProbNNk)
{
    if((p_ProbNNp - p_ProbNNk)>0)
        return true;
    else 
        return false;
}

bool FinalSelection(bool Preselection, bool PassBDT, bool NoDDstar)
{
    if(Preselection && PassBDT && NoDDstar)
        return true;
    else
        return false;
}
'''

func_isolation = '''
double GetMuCharge(double muID)
{
    return -muID/13;
}

double GetMass(double Lb_ISOLATION_PIDK, double Lb_ISOLATION_CHARGE, double muCharge, double PIDdiff)
{
    double m=0;
    double m_pi = 139.57018;
    double m_K = 493.677;
    if (Lb_ISOLATION_PIDK>4.&&(Lb_ISOLATION_CHARGE==-muCharge || (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff<0)))
        m=m_K;
    else
        m=m_pi;
    return m;
}

bool isKenriched(double Lb_ISOLATION_BDT, double ISOBDTcut, double Lb_ISOLATION_BDT2, double ISOBDT2cut, double mLc12, double mTOT, double m1, double m2)
{
    double m_K = 493.677;
    if(Lb_ISOLATION_BDT>ISOBDTcut && Lb_ISOLATION_BDT2>ISOBDT2cut)
    {
        if (mLc12>2700 && mTOT<5620 && (m1==m_K || m2==m_K))
            return true;
        else
            return false;
    }
    else
        return false;
}
'''

def CheckIsolation(df,inputFile,ISOBDTcut, ISOBDT2cut):
    #Masses pi, K, p, mu, Lc
    m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
    m_K = 493.677 #+/- 0.016 MeV (PDG)
    m_p = 938.272081 #+/- 0.000006 MeV (PDG)
    m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
    m_Lc = 2286.46 #+/- 0.14 MeV (PDG)
    df1 = df.Define("E1pi","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX,2)+TMath::Power(Lb_ISOLATION_PY,2)+TMath::Power(Lb_ISOLATION_PZ,2)+"+str(m_pi)+"*"+str(m_pi)+")")
    df1 = df1.Define("E2pi","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX2,2)+TMath::Power(Lb_ISOLATION_PY2,2)+TMath::Power(Lb_ISOLATION_PZ2,2)+"+str(m_pi)+"*"+str(m_pi)+")")
    df1 = df1.Define("ELc","TMath::Sqrt(TMath::Power(Lc_PX,2)+TMath::Power(Lc_PY,2)+TMath::Power(Lc_PZ,2)+"+str(m_Lc)+"*"+str(m_Lc)+")")
    df1 = df1.Define("pLc12_x","Lc_PX+Lb_ISOLATION_PX+Lb_ISOLATION_PX2")
    df1 = df1.Define("pLc12_y","Lc_PY+Lb_ISOLATION_PY+Lb_ISOLATION_PY2")
    df1 = df1.Define("pLc12_z","Lc_PZ+Lb_ISOLATION_PZ+Lb_ISOLATION_PZ2")
    df1 = df1.Define("mLc12","TMath::Sqrt(TMath::Power(ELc+E1pi+E2pi,2)+TMath::Power(pLc12_x,2)+TMath::Power(pLc12_y,2)+TMath::Power(pLc12_z,2))")
    df1 = df1.Define("muCharge","GetMuCharge(mu_ID)")
    df1 = df1.Define("PIDdiff","Lb_ISOLATION_PIDp - Lb_ISOLATION_PIDK")
    df1 = df1.Define("PIDdiff2","Lb_ISOLATION_PIDp2 - Lb_ISOLATION_PIDK2")
    df1 = df1.Define("m1","GetMass(Lb_ISOLATION_PIDK, Lb_ISOLATION_CHARGE, muCharge,PIDdiff)")
    df1 = df1.Define("m2","GetMass(Lb_ISOLATION_PIDK2, Lb_ISOLATION_CHARGE2, muCharge,PIDdiff2)")
    df1 = df1.Define("E1","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX,2)+TMath::Power(Lb_ISOLATION_PY,2)+TMath::Power(Lb_ISOLATION_PZ,2)+TMath::Power(m1,2))")
    df1 = df1.Define("E2","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX2,2)+TMath::Power(Lb_ISOLATION_PY2,2)+TMath::Power(Lb_ISOLATION_PZ2,2)+TMath::Power(m2,2))")
    df1 = df1.Define("Emu","TMath::Sqrt(TMath::Power(mu_PX,2)+TMath::Power(mu_PY,2)+TMath::Power(mu_PZ,2)+"+str(m_mu)+"*"+str(m_mu)+")")
    df1 = df1.Define("Etot","ELc+E1+E2+Emu")
    df1 = df1.Define("pTOT_x","pLc12_x+mu_PX")
    df1 = df1.Define("pTOT_y","pLc12_y+mu_PY")
    df1 = df1.Define("pTOT_z","pLc12_z+mu_PZ")
    df1 = df1.Define("mTOT","TMath::Sqrt(TMath::Power(Etot,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)")
    df1 = df1.Define("isKenriched","isKenriched(Lb_ISOLATION_BDT,"+str(ISOBDTcut)+", Lb_ISOLATION_BDT,"+str(ISOBDT2cut)+", mLc12, mTOT, m1,  m2)")
    df1 = df1.Define("isIsolated","Lb_ISOLATION_BDT<"+str(ISOBDTcut))
    column_name_vector = r.std.vector('string')()
    column_name_vector.push_back("isKenriched")
    column_name_vector.push_back("isIsolated")
    df1.Snapshot("DecayTree",inputFile[0:-5]+'_Isolation.root',column_name_vector)
    return


if __name__ == "__main__":
     ISOBDTcut =0.35
     ISOBDT2cut=0.2
     
     r.gInterpreter.Declare(func_code)
     r.gInterpreter.Declare(func_isolation)

     dtype = ['Data','DataSS','FakeMu','FakeMuSS']
     polarities=['MagUp','MagDown']
     for polarity in polarities:
        print(polarity)
        for dt in dtype:
            print('   ', dt)
            inputFile = filedir+'Lb_'+dt+'_'+polarity+'.root'
            f=r.TFile(inputFile)
            tname = 'tupleout/DecayTree'
            t=f.Get(tname)

            tname2 = 'DecayTree'
            bdtfname = GetBDTfile(dt,polarity)
            fbdt = r.TFile(bdtfname)
            tbdt = fbdt.Get(tname2)
            t.AddFriend(tbdt)

            df0 = r.RDataFrame(t)

            CheckIsolation(df0,inputFile,ISOBDTcut, ISOBDT2cut)
            fIso = r.TFile(inputFile[0:-5]+'_Isolation.root')
            tIso = fIso.Get(tname2)
            t.AddFriend(tIso)

            df0 = r.RDataFrame(t)
            df1 = df0.Define("L0","L0Trigger(Lb_L0Global_TIS,Lb_L0HadronDecision_TOS)")
            df1 = df1.Define("HLT1","HLT1Trigger(Lc_Hlt1TrackMVADecision_TOS,Lc_Hlt1TwoTrackMVADecision_TOS)")

            if dt=='Data' or dt=='DataSS':
                df1 = df1.Define("HLT2","HLT2Trigger(Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS)")
            else:
                df1 = df1.Define("HLT2","HLT2Trigger(Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS)")

            df1 = df1.Define("Trigger","PassTrigger(L0, HLT1, HLT2)")
            df2 = df1.Define("LcMass","LcMassCut(Lc_M)")
            df2 = df2.Define("PIDCalib", "ApplyPIDCalibCuts(nTracks, K_P, K_PT, p_P, p_PT,pi_P, pi_PT, mu_P, mu_PT)")
            if dt=='Data' or dt=='DataSS':
                df2 = df2.Define("MuCuts","ApplyMuCuts(mu_PIDmu, mu_PIDK, mu_PIDp, mu_PIDe)")   
            else: 
                df2 = df2.Define("MuCuts","true")
            df2 = df2.Define("Preselection","GetFinalPreselection(Trigger, PIDCalib, LcMass, MuCuts)")
            df3 = df2.Define("PassBDT","PassBDT(bdt)")
            df3 = df3.Define("NoDDstar","RemoveDDstar(p_ProbNNp, p_ProbNNk)")
            df3 = df3.Define("FinalSel","FinalSelection(Preselection, PassBDT,NoDDstar)")
            
            column_name_vector = r.std.vector('string')()
            column_name_vector.push_back("L0")
            column_name_vector.push_back("HLT1")
            column_name_vector.push_back("HLT2")
            column_name_vector.push_back("Trigger")
            column_name_vector.push_back("LcMass")
            column_name_vector.push_back("PIDCalib")
            column_name_vector.push_back("MuCuts")
            column_name_vector.push_back("Preselection")
            column_name_vector.push_back("PassBDT")
            column_name_vector.push_back("NoDDstar")
            column_name_vector.push_back("FinalSel")
            column_name_vector.push_back("isKenriched")
            column_name_vector.push_back("isIsolated")

            print("Writing the preselection tree ...")
            df3.Snapshot("DecayTree",inputFile[0:-5]+'_preselectionVars.root',column_name_vector)
            print("Tree written.")
            
            df4 = df3.Filter("FinalSel==true&&isIsolated==true")
            print("Writing the Isolated tree ...")
            df4.Snapshot("DecayTree",inputFile[0:-5]+'_preselected_iso.root')
            print("Tree written.")
            print("Adding Sweights")
            ofname = ComputeSweights(inputFile[0:-5]+'_preselected_iso.root',dt,polarity,'iso')
            print('Created: '+ofname)

            df5 = df3.Filter("FinalSel==true&&isKenriched==true")
            print("Writing the Kenriched tree ...")
            df5.Snapshot("DecayTree",inputFile[0:-5]+'_preselected_Kenr.root')
            print("Tree written.")
            ofname = ComputeSweights(inputFile[0:-5]+'_preselected_Kenr.root',dt,polarity,'Kenriched') 
            print('Created: '+ofname)




