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
    bdtfname = filedir+'Lb_'+dtype+'_'+polarity+'_MVA.root'
    if os.path.isfile(bdtfname):
        print('BDT file already created')
    else:
        print()
        print('>>>   Creating file with BDT variable')
        print()
        AddBDTinfo(ifname, 'DecayTree', bdtfname, 'MC',
                   pickled_model_path = '../PIDGen_PIDCalib_MVA/xgb_reg.pkl')
    return bdtfname


r.gInterpreter.Declare(func_code)
inputFile = filedir+'Lb_Lcmunu_MagUp_PID.root'

f=r.TFile(inputFile)
tname = 'DecayTree'
t=f.Get(tname)
print(t)


bdtfname = GetBDTfile('Lcmunu','MagUp')
fbdt = r.TFile(bdtfname)
tbdt = fbdt.Get(tname)
t.AddFriend(tbdt)

df0 = r.RDataFrame(t)
df1 = df0.Define("TruthMatch","TruthMatchLb(Lc_BKGCAT, Lb_BKGCAT, Lb_TRUEID, Lc_TRUEID, mu_TRUEID, mu_MC_MOTHER_ID)")
df1 = df1.Define("L0True","L0Trigger(Lc_L0HadronDecision_TOS,Lb_L0Global_TIS)")
df1 = df1.Define("HLT1True","HLT1Trigger(Lc_Hlt1TrackMVADecision_TOS, Lc_Hlt1TwoTrackMVADecision_TOS)")
df1 = df1.Define("TriggerTrue","PassTrigger(L0True, HLT1True,Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS)")
df2 = df1.Define("LcMass","LcMassCut(Lc_M)")
df2 = df2.Define("PIDCalib", "ApplyPIDCalibCuts(nTracks, K_P, K_PT, p_P, p_PT,pi_P, pi_PT, mu_P, mu_PT)")
df2 = df2.Define("PreselectionTrue","GetFinalPreselection(TruthMatch,TriggerTrue,PIDCalib, LcMass)")
df3 = df2.Define("PassBDT","PassBDT(bdt)")
df3 = df3.Define("FinalSelTrue","FinalSelection(PreselectionTrue, PassBDT)")


filtered = []
filtered.append(df3.Filter('TruthMatch==true', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true', 'CutL0_true'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true','CutL0&HLT1_True'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS','TriggerTrue'))
filtered.append(df3.Filter('TruthMatch==true&&TriggerTrue==true&&LcMass==true','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&TriggerTrue==true&&LcMass==true&&PIDCalib==true','CutPresel'))
filtered.append(df3.Filter('PreselectionTrue==true&&PassBDT==true','CutBDT'))
print()
print('Using TRUE Trigger variables')
for filt in filtered:
    rep = filt.Report()
    rep.Print()


#Adding Lc_PT>2000 cut
filtered = []
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true', 'CutL0_true'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&HLT1True==true','CutL0&HLT1_True'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS&&Lc_PT>2000','TriggerTrue'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&PIDCalib==true','CutPresel'))
filtered.append(df3.Filter('Lc_PT>2000&&PreselectionTrue==true&&PassBDT==true','CutBDT'))
print()
print('Using TRUE Trigger variables AND Lc_PT>2000 cut')
for filt in filtered:
    rep = filt.Report()
    rep.Print()

#Adding PIDcuts
filtered = []
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0', 'CutL0_true'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&HLT1True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0','CutL0&HLT1_True'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0','TriggerTrue'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&PIDCalib==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0','CutPresel'))
filtered.append(df3.Filter('Lc_PT>2000&&PreselectionTrue==true&&PassBDT==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0','CutBDT'))
print()
print('Using TRUE Trigger variables AND Lc_PT>2000 cut AND PID cuts')
for filt in filtered:
    rep = filt.Report()
    rep.Print()

#Adding Flight Lenght CHI2 cut
filtered = []
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.', 'CutL0_true'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&HLT1True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.','CutL0&HLT1_True'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.','TriggerTrue'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&PIDCalib==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.','CutPresel'))
filtered.append(df3.Filter('Lc_PT>2000&&PreselectionTrue==true&&PassBDT==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.','CutBDT'))
print()
print('Using TRUE Trigger variables AND Lc_PT>2000 cut AND PID cuts AND FDCHI2>50 cut')
for filt in filtered:
    rep = filt.Report()
    rep.Print()

#Adding daughterPcut >5 GeV
filtered = []
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000', 'CutL0_true'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&HLT1True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000','CutL0&HLT1_True'))
filtered.append(df3.Filter('TruthMatch==true&&L0True==true&&HLT1True==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS&&Lc_PT>2000&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000','TriggerTrue'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&Lc_PT>2000&&TriggerTrue==true&&LcMass==true&&PIDCalib==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000','CutPresel'))
filtered.append(df3.Filter('Lc_PT>2000&&PreselectionTrue==true&&PassBDT==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&p_P>5000&&K_P>5000&&pi_P>5000','CutBDT'))
print()
print('Using TRUE Trigger variables AND Lc_PT>2000 cut AND PID cuts AND FDCHI2>50 cut AND Daughter momentum > 5000 MeV')
for filt in filtered:
    rep = filt.Report()
    rep.Print()

# Create File for Patrick
'''
df4 = df3.Filter('TruthMatch==true&&Lc_PT>2000&&L0True==true&&HLT1True==true&&K_PIDK>2&&pi_PIDK<4&&p_PIDp>0&&Lb_FDCHI2_OWNPV>50.&&LcMass==true&&PIDCalib==true&&PassBDT==true')
print("Writing New tree ...")
df4.Snapshot("DecayTree",inputFile[0:-5]+'_forPatrick.root')
print("Tree written.")
'''
