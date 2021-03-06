import uproot
import sys, os
import pandas as pd
import numpy as np
from root_pandas import to_root
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo
from AddSweights import *
from argparse import ArgumentParser

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'

def EmulateTrigger(inputFile,restartL0HTOS,restartL0GTIS,restartHLT1_1,restartHLT1_2):
    os.chdir('../TrackerOnlyEmulators/')
    if restartL0HTOS: 
        rL0TOS = '--L0HTOS'
    else:
        rL0TOS = ''
    if restartL0GTIS: 
        rL0TIS = '--L0GTIS'
    else:
        rL0TIS = ''
    if restartHLT1_1: 
        rHLT11 = '--HLT1_1'
    else:
        rHLT11 = ''
    if restartHLT1_2: 
        rHLT12 = '--HLT1_2'
    else:
        rHLT12 = ''
    os.system('python run_triggers_emulators.py -f %s %s %s %s %s' %(inputFile,rL0TOS,rL0TIS,rHLT11,rHLT12))
    os.chdir('../Preselection/')
    return

def GetBDTfile(dtype,polarity):
    ifname = filedir+'Lb_'+dtype+'_'+polarity+'_full.root'
    bdtfname = ifname[0:-5]+'_MVA.root'
    if os.path.isfile(bdtfname):
        print('BDT file already created')
    else:
        print()
        print('>>>   Creating file with BDT variable')
        print()
        AddBDTinfo(ifname, 'tupleout/DecayTree', bdtfname, 'MC',
                   pickled_model_path = '../PIDGen_PIDCalib_MVA/xgb_reg.pkl')
    return bdtfname


func_code = '''
bool TruthMatchCharm(double Lc_BKGCAT, double Lb_BKGCAT, double Lb_TRUEID, double Lc_TRUEID, double mu_TRUEID, double mu_MC_MOTHER_ID)
{
    if(Lc_BKGCAT<30 && Lb_BKGCAT<50)
    {
        if(abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13)
        {
            if(abs(mu_MC_MOTHER_ID)==431 || abs(mu_MC_MOTHER_ID)==421 || abs(mu_MC_MOTHER_ID)==411)
            return true;
        }
    }
    return false;
}

bool TruthMatchTau(double Lc_BKGCAT, double Lb_BKGCAT, double Lb_TRUEID, double Lc_TRUEID, double mu_TRUEID, double mu_MC_MOTHER_ID)
{
    if(Lc_BKGCAT<30 && Lb_BKGCAT<50)
    {
        if(abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13)
        {
            if(abs(mu_MC_MOTHER_ID)==15)
                return true;
        }
    }
    return false;
}

bool TruthMatchLb(double Lc_BKGCAT, double Lb_BKGCAT, double Lb_TRUEID, double Lc_TRUEID, double mu_TRUEID, double mu_MC_MOTHER_ID)
{
    if(Lc_BKGCAT<30 && Lb_BKGCAT<50)
    {
        if(abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13)
        {
            if(abs(mu_MC_MOTHER_ID)==5122)
                return true;
        }
    }
    return false;
}

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

bool PassTrigger(bool L0, bool HLT1)
{
    //cout<< "L0: "<<L0<< "  HLT1:" << HLT1<< " HLT2: "<<HLT2<<endl;
    if(L0 && HLT1)
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

bool GetFinalPreselection(bool TruthMatch, bool Trigger, bool PIDCalibCuts, bool LcMass)
{
    if (TruthMatch && Trigger && PIDCalibCuts && LcMass)
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

bool FinalSelection(bool Preselection, bool PassBDT)
{
    if(Preselection && PassBDT)
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

def PlotFitTemplate(df,dtype,polarity):
    h_El = df.Histo1D(("h_El", "E*_{#mu}", 10, 0., 2600.)
,"FitVar_El_mLc","Event_PIDCalibEffWeight")
    h_q2 = df.Histo1D(("h_q2", "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","Event_PIDCalibEffWeight");
    h_Mmiss2 = df.Histo1D(("h_Mmiss2", "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","Event_PIDCalibEffWeight");
    c = r.TCanvas("c","c",1500,500)
    c.Divide(3,1)
    c.cd(1)
    h_El.Draw("hist")
    c.cd(2)
    h_q2.Draw("hist")
    c.cd(3)
    h_Mmiss2.Draw("hist")
    c.SaveAs("plots/FullMC/"+dtype+'_'+polarity+'_FitVars.png')
    return 
    

def PrintCutsEfficiencies(df):
    filtered = []
    filtered.append(df.Filter('TruthMatch==true', 'TruthMatch'))
    filtered.append(df.Filter('TruthMatch==true&&L0==true', 'CutL0'))
    filtered.append(df.Filter('TruthMatch==true&&HLT1==true', 'CutHLT1'))
    filtered.append(df.Filter('TruthMatch==true&&Trigger==true','CutL0&HLT1'))
    filtered.append(df.Filter('TruthMatch==true&&Trigger==true&&LcMass==true','CutLcM'))
    filtered.append(df.Filter('TruthMatch==true&&Trigger==true&&LcMass==true&&PIDCalib==true','CutPresel'))
    #filtered.append(df.Filter('PassBDT==true','CutBDT'))
    filtered.append(df.Filter('Preselection==true&&PassBDT==true','CutBDT'))
    filtered.append(df.Filter('FinalSel==true&&isIsolated==true','Isolated'))
    filtered.append(df.Filter('FinalSel==true&&isKenriched==true','Kenriched'))
    print('All stats:')
    for filt in filtered:
        r = filt.Report()
        r.Print()
    return



if __name__ == "__main__":
    ISOBDTcut =0.35
    ISOBDT2cut=0.2

    r.gInterpreter.Declare(func_code)
    r.gInterpreter.Declare(func_isolation)
    
    #dtype = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625Ds']
    dtype = ['Lcmunu']
    polarities=['MagUp','MagDown']
    #polarities=['MagUp']
    for polarity in polarities:
        print(polarity)
        for dt in dtype:
            print('   ', dt)
            inputFile = filedir+'Lb_'+dt+'_'+polarity+'_full.root'
            pidgenFile = inputFile[0:-5]+'_PIDGen.root'
            pidcalibFile = inputFile[0:-5]+'_PIDCalib.root'
            
            L0TOSFile = inputFile[0:-5]+'_wL0TOSEmulation.root'
            L0TISFile = inputFile[0:-5]+'_wL0TISEmulation.root'
            HLT1_1Trk_File = inputFile[0:-5]+'_wHLT1OneTrackEmulation_Df.root'
            HLT1_2Trk_File = inputFile[0:-5]+'_wHLT1TwoTracksEmulation.root'


            f=r.TFile(inputFile)
            tname = 'tupleout/DecayTree'
            t=f.Get(tname)

            fpidgen = r.TFile(pidgenFile)
            tname2 = 'DecayTree'
            tgen = fpidgen.Get(tname2)
            
            fpidcalib = r.TFile(pidcalibFile)
            tcalib = fpidcalib.Get(tname)

            fL0TOS = r.TFile(L0TOSFile)
            tL0TOS = fL0TOS.Get(tname2)

            fL0TIS = r.TFile(L0TISFile)
            tL0TIS = fL0TIS.Get(tname2)

            fHLT1_1 = r.TFile(HLT1_1Trk_File)
            tHLT1_1 = fHLT1_1.Get(tname2)

            fHLT1_2 = r.TFile(HLT1_2Trk_File)
            tHLT1_2 = fHLT1_2.Get(tname2)
            
            
            t.AddFriend(tcalib)
            t.AddFriend(tgen)
            t.AddFriend(tL0TOS)
            t.AddFriend(tL0TIS)
            t.AddFriend(tHLT1_1)
            t.AddFriend(tHLT1_2)
            
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

            if dt in ['Lb_LcDs','Lb_Lc2593Ds','Lb_Lc2625Ds']:
                df1 = df0.Define("TruthMatch","TruthMatchCharm(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)") 
            elif dt in ['Lb_Lctaunu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
                df1 = df0.Define("TruthMatch","TruthMatchTau(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
            else:
                df1 = df0.Define("TruthMatch","TruthMatchLb(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
            df1 = df1.Define("L0","L0Trigger(Lc_L0Hadron_TOS_emulated,Lb_L0Global_TIS_emulated)")
            df1 = df1.Define("HLT1","HLT1Trigger(Lc_HLT1TrackMVA_Emu_EffCorrected_TOS, Lc_HLT1TwoTrackMVA_Emu_TOS)")
            df1 = df1.Define("Trigger","PassTrigger(L0, HLT1)")
            df2 = df1.Define("LcMass","LcMassCut(Lc_M)")
            df2 = df2.Define("PIDCalib", "ApplyPIDCalibCuts(nTracks, K_P, K_PT, p_P, p_PT,pi_P, pi_PT, mu_P, mu_PT)")
            df2 = df2.Define("Preselection","GetFinalPreselection(TruthMatch,Trigger,PIDCalib, LcMass)")
            df3 = df2.Define("PassBDT","PassBDT(bdt)")
            df3 = df3.Define("FinalSel","FinalSelection(Preselection, PassBDT)")

            column_name_vector = r.std.vector('string')()
            column_name_vector.push_back("TruthMatch")
            column_name_vector.push_back("L0")
            column_name_vector.push_back("HLT1")
            column_name_vector.push_back("Trigger")
            column_name_vector.push_back("LcMass")
            column_name_vector.push_back("PIDCalib")
            column_name_vector.push_back("Preselection")
            column_name_vector.push_back("PassBDT")
            column_name_vector.push_back("FinalSel")
            column_name_vector.push_back("isKenriched")
            column_name_vector.push_back("isIsolated")

            print("Writing the preselection tree ...")
            df3.Snapshot("DecayTree",inputFile[0:-5]+'_preselectionVars.root',column_name_vector)
            print("Tree written.")
           
            
            PrintCutsEfficiencies(df3)
            #PlotFitTemplate(df4, dt, polarity)



