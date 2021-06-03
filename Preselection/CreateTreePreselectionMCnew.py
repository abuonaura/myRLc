'''
Author: Annarita Buonaura
Date: May 21, 2020

Description: Script to apply preselections to MC samples

How to run:
    - On full simulation (all sample, producing from scratch also the trigger files):
    python -i CreateTreePreselectionMC.py --MCfull --L0GTIS --L0HTOS --HLT1_1 --HLT1_2 all all --presel

    - On TrackerOnly sim (all sample, producing from scratch also the trigger files):
    python -i CreateTreePreselectionMC.py --MCTrackerOnly --L0GTIS --L0HTOS --HLT1_1 --HLT1_2 all all --presel

    - to see only selection efficiencies:
    python -i CreateTreePreselectionMC.py --MCTrackerOnly --L0GTIS --L0HTOS --HLT1_1 --HLT1_2 all all --eff

    - to select only 1 polarity (e.g. MagUp):
    python -i CreateTreePreselectionMC.py --MCTrackerOnly --L0GTIS --L0HTOS --HLT1_1 --HLT1_2 all MagUp --presel

    - to select only 1 sample (e.g. Lb_Lcmunu):
    python -i CreateTreePreselectionMC.py --MCTrackerOnly --L0GTIS --L0HTOS --HLT1_1 --HLT1_2 Lb_Lcmunu all --presel
'''
#import uproot
import sys, os
import pandas as pd
import numpy as np
import ROOT as r
#from root_pandas import to_root
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo
from AddSweights import *
from argparse import ArgumentParser
from LbCorr import correctLb



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

def GetBDTfile(ifname):
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
            else
                return false;
        }
        else
            return false;
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
            else
                return false;
        }
        else
            return false;
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
            else
                return false;
        }
        else
            return false;
    }
    return false;
}

bool TruthMatchBLcpbar(double Lc_BKGCAT, double Lb_BKGCAT, double Lb_TRUEID, double Lc_TRUEID, double mu_TRUEID, double mu_MC_MOTHER_ID)
{
    if(Lc_BKGCAT<30 && Lb_BKGCAT<50)
    {
        if(abs(Lb_TRUEID)==521&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13)
        {
            if(abs(mu_MC_MOTHER_ID)==521)
                return true;
            else
                return false;
        }
        else
            return false;
    }
    return false;
}

bool PassHLT2Selections(double Lc_PT, double K_PIDK, double pi_PIDK, double p_PIDp, double Lb_FDCHI2_OWNPV, double p_P, double pi_P, double K_P)
{
    if(Lc_PT>2000&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000)
        return true;
    return false;
}

bool L0Trigger(bool L0TOS, bool L0TIS, int nSPDHits)
{
    if (nSPDHits<450)
    {
        if(L0TOS || L0TIS)
            return true;
        else
            return false;
    }
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

bool GetFinalPreselection(bool TruthMatch, bool MoreHLT2Sel, bool Trigger, bool PIDCalibCuts, bool LcMass)
{
    if (TruthMatch && MoreHLT2Sel && Trigger && PIDCalibCuts && LcMass)
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
    double m_pi = 139.57018;
    double m_K = 493.677;
    double m_p = 938.27208;
    double m = m_pi;
    if (Lb_ISOLATION_PIDK>4.)
    {
        if (Lb_ISOLATION_CHARGE==-muCharge || (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff<=0.))
            m=m_K;
        else if (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff>0)
            m=m_p;
    }
    return m;
}

bool isKenriched(double Lb_ISOLATION_BDT, double ISOBDTcut, double Lb_ISOLATION_BDT2, double ISOBDT2cut, double mLc12, double mTOT, double m1, double m2, double Type, double Type2, double NNghost, double NNghost2)
{
    double m_K = 493.677;
    if (Type==3 && NNghost<0.2 && Type2==3 && NNghost2<0.2)
    {
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
    else
      return false;
}

bool isLcpipi(bool isKenriched, double Lb_ISOLATION_BDT, double ISOBDTcut, double Lb_ISOLATION_BDT2, double ISOBDT2cut, double mLc12, double mTOT, double Type, double Type2, double NNghost, double NNghost2, double charge, double charge2)
{
    if (Type==3 && NNghost<0.2 && Type2==3 && NNghost2<0.2)
    {
        if(Lb_ISOLATION_BDT>ISOBDTcut && Lb_ISOLATION_BDT2>ISOBDT2cut && isKenriched==false && charge==-charge2)
        {
            if(mLc12<2700)
            {
                return true;
                }
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}
'''

func_doublecharm = '''
bool SelectMBody(double Lb_TrueHadron_D0_ID, double Lb_TrueHadron_D1_ID, double Lb_TrueHadron_D2_ID, double D0pdg)
{
    if (abs(Lb_TrueHadron_D0_ID)==D0pdg && Lb_TrueHadron_D1_ID!=0 && Lb_TrueHadron_D2_ID!=0)
        return true;
    else
        return false;
}

bool Select2Body(double Lb_TrueHadron_D0_ID, double Lb_TrueHadron_D1_ID, double Lb_TrueHadron_D2_ID, double D0pdg)
{
    if (abs(Lb_TrueHadron_D0_ID)==D0pdg && Lb_TrueHadron_D1_ID!=0 && Lb_TrueHadron_D2_ID==0)
        return true;
    else
        return false;
}

double ComputeMbodyWeight(bool mbody, double Lb_TrueHadron_D0_PX, double Lb_TrueHadron_D0_PY, double Lb_TrueHadron_D0_PZ, double Lb_TrueHadron_D0_PE, double Lb_TrueHadron_D1_PX, double Lb_TrueHadron_D1_PY, double Lb_TrueHadron_D1_PZ, double Lb_TrueHadron_D1_PE)
{
    if (mbody==true)
    {
        TLorentzVector Lc_LV(Lb_TrueHadron_D0_PX,Lb_TrueHadron_D0_PY,Lb_TrueHadron_D0_PZ, Lb_TrueHadron_D0_PE);
        TLorentzVector D_LV (Lb_TrueHadron_D1_PX,Lb_TrueHadron_D1_PY,Lb_TrueHadron_D1_PZ, Lb_TrueHadron_D1_PE);
        double mLcD2 = (Lc_LV+D_LV).M2();
        double mLc = Lc_LV.M();
        double mD = D_LV.M();
        double mLb = 5619.6;
        double mK = 493.677;
        double weight = (mLcD2 - TMath::Power(mLc+mD,2))/(TMath::Power(mLb-mK,2)- TMath::Power(mLc+mD,2));
        //cout<<weight<<endl;
        if (weight>0)
            weight = TMath::Sqrt(weight)-1./2;
        else
        {
            cout<<" Negative mbody weight "<<endl;
            weight = -1000.;
        }
        return weight;
    }
    else
        return 1;
}
'''


def CheckIsolation(df,inputFile,dt,ISOBDTcut, ISOBDT2cut):
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
    df1 = df1.Define("mLc12","TMath::Sqrt(TMath::Power(ELc+E1pi+E2pi,2)-(TMath::Power(pLc12_x,2)+TMath::Power(pLc12_y,2)+TMath::Power(pLc12_z,2)))")
    #h = df1.Histo1D(('h_mLc12','',100,2000,5000),'mLc12')
    #h.Draw()
    #h.SaveAs('plots/mLc12_'+dt+'.png')
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
    #Lcpipi tot energy (both anti-iso particle are pions)
    df1 = df1.Define("Etot1","ELc+E1pi+E2pi+Emu")
    df1 = df1.Define("mTOT1","TMath::Sqrt(TMath::Power(Etot1,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)")
    df1 = df1.Define("isKenriched","isKenriched(Lb_ISOLATION_BDT,"+str(ISOBDTcut)+", Lb_ISOLATION_BDT2,"+str(ISOBDT2cut)+", mLc12, mTOT, m1,  m2, Lb_ISOLATION_Type, Lb_ISOLATION_Type2, Lb_ISOLATION_NNghost, Lb_ISOLATION_NNghost2)")
    df1 = df1.Define("isLcpipi","isLcpipi(isKenriched,Lb_ISOLATION_BDT,"+str(ISOBDTcut)+", Lb_ISOLATION_BDT2,"+str(ISOBDT2cut)+", mLc12, mTOT1, Lb_ISOLATION_Type, Lb_ISOLATION_Type2, Lb_ISOLATION_NNghost, Lb_ISOLATION_NNghost2, Lb_ISOLATION_CHARGE, Lb_ISOLATION_CHARGE2)")
    df1 = df1.Define("isIsolated","Lb_ISOLATION_BDT<"+str(ISOBDTcut))
    column_name_vector = r.std.vector('string')()
    column_name_vector.push_back("isKenriched")
    column_name_vector.push_back("isIsolated")
    column_name_vector.push_back("isLcpipi")
    df1.Snapshot("DecayTree",inputFile[0:-5]+'_Isolation.root',column_name_vector)
    return

def ApplyLbCorrection(df,inputFile,dt):
    df1 = df.Define("w_LbCorr","GetLbCorrWeight(Lb_TRUEP_X, Lb_TRUEP_Y, Lb_TRUEP_Z, Lb_TRUEPT)")
    column_name_vector = r.std.vector('string')()
    column_name_vector.push_back("w_LbCorr")
    df1.Snapshot("DecayTree",inputFile[0:-5]+'_LbCorrection.root',column_name_vector)
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
    filtered.append(df.Filter('TruthMatch==true&&MoreHLT2Sel==true', 'AdditionalHLT2Selections'))
    filtered.append(df.Filter('TruthMatch==true&&MoreHLT2Sel==true&&L0==true', 'CutL0'))
    filtered.append(df.Filter('TruthMatch==true&&MoreHLT2Sel==true&&L0==true&&HLT1==true', 'CutHLT1'))
    filtered.append(df.Filter('TruthMatch==true&&MoreHLT2Sel==true&&Trigger==true&&LcMass==true','CutLcM'))
    filtered.append(df.Filter('TruthMatch==true&&MoreHLT2Sel==true&&Trigger==true&&LcMass==true&&PIDCalib==true','CutPIDCalib->Preselection'))
    filtered.append(df.Filter('Preselection==true&&PassBDT==true','FinalCut'))
    print('All stats:')
    for filt in filtered:
        r = filt.Report()
        r.Print()
    return



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    parser.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    parser.add_argument('--L0HTOS',dest='restartL0HTOS', help="Forces to reproduce L0Hadron TOS file", required=False, default=False, action='store_true')
    parser.add_argument('--L0GTIS',dest='restartL0GTIS', help="Forces to reproduce L0Global TIS file", required=False, default=False, action='store_true')
    parser.add_argument('--HLT1_1',dest='restartHLT1_1', help="Forces to reproduce HLT1 one track file", required=False, default=False, action='store_true')
    parser.add_argument('--HLT1_2',dest='restartHLT1_2', help="Forces to reproduce HLT1 two track file", required=False, default=False, action='store_true')
    parser.add_argument('datatype', help = 'which mc sample we want to run on', default = 'all')
    parser.add_argument('polarity',choices=['MagUp','MagDown','all'], help = 'which data sample we want to run on', default = 'all')
    parser.add_argument('--presel',dest='preselect',help='If this option is inserted the preselection file + the iso/Kenriched sample files are produced', default=False, action='store_true')
    parser.add_argument('--eff',dest='efficiency',help='If this option is inserted, only the efficiency are printed out', default=False, action='store_true')

    options = parser.parse_args()
    MCfull = options.MCfull
    MCTO = options.MCTO
    restartL0HTOS = options.restartL0HTOS
    restartL0GTIS = options.restartL0GTIS
    restartHLT1_1 = options.restartHLT1_1
    restartHLT1_2 = options.restartHLT1_2
    datatype=options.datatype
    polarity=options.polarity
    preselect = options.preselect
    efficiencies = options.efficiency
    
    if MCfull==True:
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_trueTrigger/'
        #filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
    if MCTO==True:
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/'

    ISOBDTcut =0.35
    ISOBDT2cut=0.2

    r.gInterpreter.Declare(func_code)
    r.gInterpreter.Declare(func_isolation)
    r.gInterpreter.Declare(func_doublecharm)
    
    if(MCfull):
        datatypes = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds']
    if(MCTO):
        datatypes = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','B_Lcpbarmunu','Lb_Lc2765munu','Lb_Lc2880munu']
    polarities=['MagUp','MagDown']

    if datatype!='all':
        datatypes=[datatype]
    if polarity!='all':
        polarities=[polarity]


    for polarity in polarities:
        print(polarity)
        for dt in datatypes:
            print('   ', dt)
            if MCfull==True:
                inputFile = filedir+dt+'_'+polarity+'_full.root'
            if MCTO==True:
                inputFile = filedir+dt+'_'+polarity+'.root'
            pidgenFile = inputFile[0:-5]+'_PIDGen.root'
            pidcalibFile = inputFile[0:-5]+'_PIDCalib.root'
            EmulateTrigger(inputFile,restartL0HTOS,restartL0GTIS,restartHLT1_1,restartHLT1_2)

            L0TOSFile = inputFile[0:-5]+'_wL0TOSEmulation.root'
            L0TISFile = inputFile[0:-5]+'_wL0TISEmulation.root'
            HLT1_File = inputFile[0:-5]+'_wHLT1Emulation.root'
            #HLT1_1Trk_File = inputFile[0:-5]+'_wHLT1OneTrackEmulation_Df.root'
            #HLT1_2Trk_File = inputFile[0:-5]+'_wHLT1TwoTracksEmulation.root'
            FFcorr_File = inputFile[0:-5]+'_FFcorrections.root'

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
            
            fHLT1 = r.TFile(HLT1_File)
            tHLT1 = fHLT1.Get(tname2)

            #fHLT1_1 = r.TFile(HLT1_1Trk_File)
            #tHLT1_1 = fHLT1_1.Get(tname2)

            #fHLT1_2 = r.TFile(HLT1_2Trk_File)
            #tHLT1_2 = fHLT1_2.Get(tname2)
            
            #add FFcorrection file only for mu/tau samples
            if dt in ['Lb_Lcmunu','Lb_Lctaunu','Lb_Lc2593munu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
                fFFcorr = r.TFile(FFcorr_File)
                tFFcorr = fFFcorr.Get(tname2)
                t.AddFriend(tFFcorr)

            
            t.AddFriend(tcalib)
            t.AddFriend(tgen)
            t.AddFriend(tL0TOS)
            t.AddFriend(tL0TIS)
            t.AddFriend(tHLT1)
            #t.AddFriend(tHLT1_1)
            #t.AddFriend(tHLT1_2)
            
            bdtfname = GetBDTfile(inputFile)
            fbdt = r.TFile(bdtfname)
            tbdt = fbdt.Get(tname2)
            t.AddFriend(tbdt)

            df0 = r.RDataFrame(t)
            

            print('Writing isolation variables')
            CheckIsolation(df0,inputFile,dt,ISOBDTcut, ISOBDT2cut)
            fIso = r.TFile(inputFile[0:-5]+'_Isolation.root')
            tIso = fIso.Get(tname2)
            t.AddFriend(tIso)

            #ApplyLbCorrection(df0,inputFile,dt)
            print('Computing Lb Corrections')
            correctLb(inputFile,tname)
            fLbCorr = r.TFile(inputFile[0:-5]+'_LbCorrection.root')
            tLbCorr = fLbCorr.Get(tname2)
            t.AddFriend(tLbCorr)

            df0 = r.RDataFrame(t)

            if dt in ['Lb_LcDs','Lb_Lc2593Ds','Lb_Lc2625Ds']:
                df1 = df0.Define("TruthMatch","TruthMatchCharm(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)") 
            if dt in ['Lb_Lctaunu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
                df1 = df0.Define("TruthMatch","TruthMatchTau(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
            if dt in ['Lb_Lcmunu','Lb_Lc2593munu','Lb_Lc2625munu','Lb_Lc2880munu','Lb_Lc2765munu']:
                df1 = df0.Define("TruthMatch","TruthMatchLb(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
            if dt=='B_Lcpbarmunu':
                df1 = df0.Define("TruthMatch","TruthMatchBLcpbar(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")

            df1 = df1.Define("MoreHLT2Sel","PassHLT2Selections(Lc_PT, K_PIDK, pi_PIDK, p_PIDp, Lb_FDCHI2_OWNPV, p_P, pi_P, K_P)")
            df1 = df1.Define("L0","L0Trigger(Lc_L0Hadron_TOS_emulated,Lb_L0Global_TIS_emulated, nSPDHits)")
            df1 = df1.Define("HLT1","HLT1Trigger(Lc_HLT1TrackMVA_Emu_EffCorrected_TOS, Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS)")
            df1 = df1.Define("Trigger","PassTrigger(L0, HLT1)")
            df2 = df1.Define("LcMass","LcMassCut(Lc_M)")
            df2 = df2.Define("PIDCalib", "ApplyPIDCalibCuts(nTracks, K_P, K_PT, p_P, p_PT,pi_P, pi_PT, mu_P, mu_PT)")
            df2 = df2.Define("Preselection","GetFinalPreselection(TruthMatch,MoreHLT2Sel,Trigger,PIDCalib, LcMass)")
            df3 = df2.Define("PassBDT","PassBDT(bdt)")
            df3 = df3.Define("FinalSel","FinalSelection(Preselection, PassBDT)")
            if dt in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
                if dt=='Lb_LcDs':
                    df3 = df3.Define("D0pdg", "4122")
                if dt=='Lb_Lc2593Ds':
                    df3 = df3.Define("D0pdg", "14122")
                if dt=='Lb_Lc2625Ds':
                    df3 = df3.Define("D0pdg", "104124")
                df3 = df3.Define("twobody","Select2Body(Lb_TrueHadron_D0_ID, Lb_TrueHadron_D1_ID, Lb_TrueHadron_D2_ID, D0pdg)")
                df3 = df3.Define("mbody","SelectMBody(Lb_TrueHadron_D0_ID, Lb_TrueHadron_D1_ID, Lb_TrueHadron_D2_ID, D0pdg)")
                df3 = df3.Define("w_mbody","ComputeMbodyWeight(mbody, Lb_TrueHadron_D0_PX, Lb_TrueHadron_D0_PY, Lb_TrueHadron_D0_PZ, Lb_TrueHadron_D0_PE, Lb_TrueHadron_D1_PX, Lb_TrueHadron_D1_PY, Lb_TrueHadron_D1_PZ, Lb_TrueHadron_D1_PE)")

            column_name_vector = r.std.vector('string')()
            column_name_vector.push_back("TruthMatch")
            column_name_vector.push_back("MoreHLT2Sel")
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
            column_name_vector.push_back("isLcpipi")
            column_name_vector.push_back("w_LbCorr")
            if dt in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
                column_name_vector.push_back("twobody")
                column_name_vector.push_back("mbody")
                column_name_vector.push_back("w_mbody")
            if dt in ['Lb_Lcmunu','Lb_Lctaunu','Lb_Lc2593munu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2625taunu']:
                column_name_vector.push_back("Event_FFcorr")

            if preselect==True:
                print('ok') 
                print("Writing the preselection tree ...")
                df3.Snapshot("DecayTree",inputFile[0:-5]+'_preselectionVars.root',column_name_vector)
                print("Tree written.")
                
                df5 = df3.Filter("FinalSel==true&&isIsolated==true")
                print("Writing the Isolated tree ...")
                df5.Snapshot("DecayTree",inputFile[0:-5]+'_preselected_iso.root')
                print("Tree written.")

                '''
                df6 = df3.Filter("FinalSel==true&&isKenriched==true")
                print("Writing the Kenriched tree ...")
                df6.Snapshot("DecayTree",inputFile[0:-5]+'_preselected_Kenr.root')
                print("Tree written.")
                
                df7 = df3.Filter("FinalSel==true&&isLcpipi==true")
                print("Writing the Lcpipi tree ...")
                df7.Snapshot("DecayTree",inputFile[0:-5]+'_preselected_Lcpipi.root')
                print("Tree written.")
                '''


            if efficiencies==True:
                PrintCutsEfficiencies(df3)
            #PlotFitTemplate(df4, dt, polarity)



