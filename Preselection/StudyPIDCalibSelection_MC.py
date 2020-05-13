import sys, os
import pandas as pd
import numpy as np
import ROOT as r
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo


filedir = '/disk/lhcb_data2/RLcMuonic2016/MC/'


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

bool PassTrigger(bool L0, bool HLT1, double HLT2)
{
    //cout<< "L0: "<<L0<< "  HLT1:" << HLT1<< " HLT2: "<<HLT2<<endl;
    if(L0 && HLT1&&HLT2)
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

def GetBDTfile(dtype,polarity):
    ifname = filedir+'Lb_'+dtype+'_'+polarity+'.root'
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


r.gInterpreter.Declare(func_code)
inputFile = filedir+'Lb_Lcmunu_MagUp.root'
L0TOSFile = inputFile[0:-5]+'_wL0TOSEmulation.root'
L0TISFile = inputFile[0:-5]+'_wL0TISEmulation.root'
HLT1_1Trk_File = inputFile[0:-5]+'_wHLT1OneTrackEmulation_Df.root'
HLT1_2Trk_File = inputFile[0:-5]+'_wHLT1TwoTracksEmulation.root'
pidgenFile = inputFile[0:-5]+'_PIDGen.root'
pidcalibFile = inputFile[0:-5]+'_PIDCalib.root'


f=r.TFile(inputFile)
tname = 'tupleout/DecayTree'
t=f.Get(tname)

tname2 = 'DecayTree'
fL0TOS = r.TFile(L0TOSFile)
tL0TOS = fL0TOS.Get(tname2)

fL0TIS = r.TFile(L0TISFile)
tL0TIS = fL0TIS.Get(tname2)

fHLT1_1 = r.TFile(HLT1_1Trk_File)
tHLT1_1 = fHLT1_1.Get(tname2)

fHLT1_2 = r.TFile(HLT1_2Trk_File)
tHLT1_2 = fHLT1_2.Get(tname2)

fpidgen = r.TFile(pidgenFile)
tname2 = 'DecayTree'
tgen = fpidgen.Get(tname2)

fpidcalib = r.TFile(pidcalibFile)
tcalib = fpidcalib.Get(tname)

t.AddFriend(tcalib)
t.AddFriend(tgen)
t.AddFriend(tL0TOS)
t.AddFriend(tL0TIS)
t.AddFriend(tHLT1_1)
t.AddFriend(tHLT1_2)

bdtfname = GetBDTfile('Lcmunu','MagUp')
fbdt = r.TFile(bdtfname)
tbdt = fbdt.Get(tname2)
t.AddFriend(tbdt)

df0 = r.RDataFrame(t)
df1 = df0.Define("TruthMatch","TruthMatchLb(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
df1 = df1.Define("L0True","L0Trigger(Lc_L0HadronDecision_TOS,Lb_L0Global_TIS)")
df1 = df1.Define("L0","L0Trigger(Lc_L0Hadron_TOS_emulated,Lb_L0Global_TIS_emulated)")
df1 = df1.Define("HLT1","HLT1Trigger(Lc_HLT1TrackMVA_Emu_EffCorrected_TOS, Lc_HLT1TwoTrackMVA_Emu_TOS)")
df1 = df1.Define("HLT1True","HLT1Trigger(Lc_Hlt1TrackMVADecision_TOS, Lc_Hlt1TwoTrackMVADecision_TOS)")
df1 = df1.Define("Trigger","PassTrigger(L0, HLT1,Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS)")
df1 = df1.Define("TriggerTrue","PassTrigger(L0True, HLT1True,Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS)")
df2 = df1.Define("LcMass","LcMassCut(Lc_M)")
df2 = df2.Define("PIDCalib", "ApplyPIDCalibCuts(nTracks, K_P, K_PT, p_P, p_PT,pi_P, pi_PT, mu_P, mu_PT)")
df2 = df2.Define("Preselection","GetFinalPreselection(TruthMatch,Trigger,PIDCalib, LcMass)")
df3 = df2.Define("PassBDT","PassBDT(bdt)")
df3 = df3.Define("FinalSel","FinalSelection(Preselection, PassBDT)")

filtered = []
filtered.append(df3.Filter('TruthMatch==true&&nTracks>0&&nTracks<700', 'nTracksCut_TM'))
filtered.append(df3.Filter('TruthMatch==true&&K_P>0&&K_P<2.E5&&K_PT>0&&K_PT<6.E4', 'KCut_TM'))
filtered.append(df3.Filter('TruthMatch==true&&pi_P>0&&pi_P<2.E5&&pi_PT>0&&pi_PT<6.E4', 'piCut_TM'))
filtered.append(df3.Filter('TruthMatch==true&&p_P>0&&p_P<2.E5&&p_PT>0&&p_PT<6.E4', 'pCut_TM'))
filtered.append(df3.Filter('TruthMatch==true&&mu_P>0&&mu_P<2.E5&&mu_PT>0&&mu_PT<6.E4', 'muCut_TM'))
print('All stats:')
for filt in filtered:
    rep = filt.Report()
    rep.Print()

df4 = df3.Filter('TruthMatch==true && Lc_PT>2000 && TriggerTrue==true && LcMass==true')
h_nTracks = df4.Histo1D('nTracks')
h_KP= df4.Histo1D('K_P')
h_KPT= df4.Histo1D('K_PT')
h_piP= df4.Histo1D('pi_P')
h_piPT= df4.Histo1D('pi_PT')
h_pP= df4.Histo1D('p_P')
h_pPT= df4.Histo1D('p_PT')
h_muP= df4.Histo1D('mu_P')
h_muPT= df4.Histo1D('mu_PT')

c = r.TCanvas('c','c')
h_nTracks.Draw()

c1 = r.TCanvas('c1','c1')
c1.Divide(2,1)
c1.cd(1)
h_KP.Draw()
c1.cd(2)
h_KPT.Draw()

c2 = r.TCanvas('c2','c2')
c2.Divide(2,1)
c2.cd(1)
h_piP.Draw()
c2.cd(2)
h_piPT.Draw()

c3 = r.TCanvas('c3','c3')
c3.Divide(2,1)
c3.cd(1)
h_pP.Draw()
c3.cd(2)
h_pPT.Draw()

c4 = r.TCanvas('c4','c4')
c4.Divide(2,1)
c4.cd(1)
h_muP.Draw()
c4.cd(2)
h_muPT.Draw()
