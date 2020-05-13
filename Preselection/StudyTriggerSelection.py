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
df2 = df2.Define("PreselectionTrue","GetFinalPreselection(TruthMatch,TriggerTrue,PIDCalib, LcMass)")
df3 = df2.Define("PassBDT","PassBDT(bdt)")
df3 = df3.Define("FinalSel","FinalSelection(Preselection, PassBDT)")
df3 = df3.Define("FinalSelTrue","FinalSelection(PreselectionTrue, PassBDT)")

filtered = []
filtered.append(df3.Filter('TruthMatch==true', 'TruthMatch'))
filtered.append(df3.Filter('TruthMatch==true&&L0==true', 'CutL0'))
filtered.append(df3.Filter('TruthMatch==true&&L0==true&&HLT1==true','CutL0&HLT1'))
filtered.append(df3.Filter('TruthMatch==true&&L0==true&&HLT1==true&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS','Trigger'))
filtered.append(df3.Filter('TruthMatch==true&&Trigger==true&&LcMass==true','CutLcM'))
filtered.append(df3.Filter('TruthMatch==true&&Trigger==true&&LcMass==true&&PIDCalib==true','CutPresel'))
filtered.append(df3.Filter('Preselection==true&&PassBDT==true','CutBDT'))

print('All stats:')
for filt in filtered:
    rep = filt.Report()
    rep.Print()


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


'''
df4 = df3.Filter("TruthMatch==true&&L0==true&&HLT1==true")
h_HLT2_0 = df4.Histo1D(("h_HLT2_0","Hlt2 after L0&&HLT1",2,0,2),"Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS")
c = r.TCanvas("c","c",500,500)
h_HLT2_0.Draw()

df5 = df3.Filter("TruthMatch==true&&L0==true&&HLT1==true&&LcMass==true&&PIDCalib==true&&PassBDT==true")
n_pass_HLT2 = df5.Filter("Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==1").Count()
n_NOTpass_HLT2 = df5.Filter("Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==0").Count()
print('Number of events after all cuts passing HLT2')
print(n_pass_HLT2.GetValue())
print('Number of events after all cuts NOT passing HLT2')
print(n_NOTpass_HLT2.GetValue())

h_HLT2_1 = df5.Histo1D(("h_HLT2_1","Hlt2 after all other cuts",2,0,2),"Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS")
c1 = r.TCanvas("c1","c1",500,500)
h_HLT2_1.Draw()

df6 = df5.Filter("Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==1")
h_El_1 = df6.Histo1D(("h_El_1", "E*_{#mu}", 10, 0., 2600.)
,"FitVar_El_mLc","Event_PIDCalibEffWeight")
h_q2_1 = df6.Histo1D(("h_q2_1", "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","Event_PIDCalibEffWeight");
h_Mmiss2_1 = df6.Histo1D(("h_Mmiss2_1", "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","Event_PIDCalibEffWeight");
c2 = r.TCanvas("c2","c2",1500,500)
c2.Divide(3,1)
c2.cd(1)
h_El_1.Draw("hist")
c2.cd(2)
h_q2_1.Draw("hist")
c2.cd(3)
h_Mmiss2_1.Draw("hist")

df7 = df5.Filter("Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==0")
h_El_0 = df7.Histo1D(("h_El_0", "E*_{#mu}", 10, 0., 2600.)
,"FitVar_El_mLc","Event_PIDCalibEffWeight")
h_q2_0 = df7.Histo1D(("h_q2_0", "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","Event_PIDCalibEffWeight");
h_Mmiss2_0 = df7.Histo1D(("h_Mmiss2_0", "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","Event_PIDCalibEffWeight");
c3 = r.TCanvas("c3","c3",1500,500)
c3.Divide(3,1)
c3.cd(1)
h_El_0.Draw("hist")
c3.cd(2)
h_q2_0.Draw("hist")
c3.cd(3)
h_Mmiss2_0.Draw("hist")
'''

'''
filedir2 = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'
inputFile = filedir2+'Lb_Lcmunu_MagUp_full.root'
preselFile = inputFile[0:-5]+'_preselectionVars.root'
pidgenFile = inputFile[0:-5]+'_PIDGen.root'
pidcalibFile = inputFile[0:-5]+'_PIDCalib.root'

tname ='tupleout/DecayTree'
f=r.TFile(inputFile)
t=f.Get(tname)

tname2 = 'DecayTree'

fpidgen = r.TFile(pidgenFile)
tname2 = 'DecayTree'
tgen = fpidgen.Get(tname2)

fpidcalib = r.TFile(pidcalibFile)
tcalib = fpidcalib.Get(tname)

fpresel = r.TFile(preselFile)
tpresel = fpresel.Get(tname2)

t.AddFriend(tcalib)
t.AddFriend(tgen)
t.AddFriend(tpresel)

n=0
nw = 0
for i in range(t.GetEntries()):
    n+=1
    t.GetEntry(i)
    tpresel.GetEntry(i)
    tcalib.GetEntry(i)
    if tpresel.FinalSel==1:
        nw+=tcalib.Event_PIDCalibEffWeight


inputFile = filedir2+'Lb_Lcmunu_MagDown_full.root'
preselFile = inputFile[0:-5]+'_preselectionVars.root'
pidgenFile = inputFile[0:-5]+'_PIDGen.root'
pidcalibFile = inputFile[0:-5]+'_PIDCalib.root'


f=r.TFile(inputFile)
t=f.Get(tname)

tname2 = 'DecayTree'


fpidgen = r.TFile(pidgenFile)
tname2 = 'DecayTree'
tgen = fpidgen.Get(tname2)

fpidcalib = r.TFile(pidcalibFile)
tcalib = fpidcalib.Get(tname)

fpresel = r.TFile(preselFile)
tpresel = fpresel.Get(tname2)

t.AddFriend(tcalib)
t.AddFriend(tgen)
t.AddFriend(tpresel)

for i in range(t.GetEntries()):
    n+=1
    t.GetEntry(i)
    tpresel.GetEntry(i)
    tcalib.GetEntry(i)
    if tpresel.FinalSel==1:
        nw+=tcalib.Event_PIDCalibEffWeight

print(n,nw)
'''
