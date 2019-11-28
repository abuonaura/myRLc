import ROOT as r
from ROOT import TMath as rm
import numpy as np
from TruthMatch import *


def GetBDTcut(BDTcut):
    bdt = r.TCut("bdt>"+str(BDTcut))
    return bdt

def GetpProbNNcut():
    pcut = r.TCut('p_ProbNNp - p_ProbNNk>0.')
    return pcut

def GetmuPIDcut():
    mucut = r.TCut('mu_PIDmu>2. && mu_PIDmu-mu_PIDK>2. && mu_PIDmu-mu_PIDp>2.')
    return mucut


def SelectMBody(tree):
    if abs(getattr(tree,'Lb_TrueHadron_D0_ID'))==4122 and getattr(tree,'Lb_TrueHadron_D1_ID')!=0 and getattr(tree,'Lb_TrueHadron_D2_ID')!=0:
        flag = True
    else:
        flag = False
    return flag

def Select2Body(tree):
    if abs(int(getattr(tree,'Lb_TrueHadron_D0_ID')))==4122 and int(getattr(tree,'Lb_TrueHadron_D1_ID'))!=0 and int(getattr(tree,'Lb_TrueHadron_D2_ID'))==0:
        flag = True
    else:
        flag = False
    return flag




def ApplyBDTcut(fname, dtype, BDTcut):
    if dtype!='MC':
        f = r.TFile(fname,'read')
        t=f.Get('DecayTree')
        tl=f.Get('LumiTuple')
        ofname = fname[0:-5]+'_BDT.root'

        of = r.TFile(ofname,'RECREATE')
        ot = r.TTree('DecayTree','DecayTree')
        otl =r.TTree('LumiTuple','LumiTuple')

        bdt = GetBDTcut(BDTcut)
        pcut = GetpProbNNcut()
        mucut = GetmuPIDcut()

        cut = r.TCut("("+bdt.GetTitle()+")")

        print('Appling cuts: ', cut)
        print('... CopyingTree with BDT cut ...')
        ot = t.CopyTree(cut.GetTitle())
        otl=tl.CloneTree()

        of.Write()
        of.Close()
        f.Close()
        print('<<<<< Created file ', ofname)
    return ofname


def ApplyFinalSelections(fname,dtype,BDTcut,mcsample):
    if dtype!='MC':
        f = r.TFile(fname,'read')
        t=f.Get('DecayTree')
        tl=f.Get('LumiTuple')
        ofname = fname[0:-5]+'_preselected.root'

        of = r.TFile(ofname,'RECREATE')
        ot = r.TTree('DecayTree','DecayTree')
        otl =r.TTree('LumiTuple','LumiTuple')

        bdt = GetBDTcut(BDTcut)
        pcut = GetpProbNNcut()
        mucut = GetmuPIDcut()

        if dtype=='Data' or dtype =='DataSS':
            cut = r.TCut("("+bdt.GetTitle()+') && ('+pcut.GetTitle()+') && ('+mucut.GetTitle()+')')
        else:
            cut = r.TCut("("+bdt.GetTitle()+') && ('+pcut.GetTitle()+')')

        print('Appling cuts: ', cut)
        print('... CopyingTree with final preselections ...')
        ot = t.CopyTree(cut.GetTitle())
        otl=tl.CloneTree()

        of.Write()
        of.Close()
        f.Close()
        print('<<<<< Created file ', ofname)
    if dtype=='MC':
        f=r.TFile(fname,'read')
        t = f.Get('DecayTree')

        ofname = fname[0:-5]+'_preselected.root'
        of = r.TFile(ofname,'RECREATE')
        ot = r.TTree('DecayTree','DecayTree')

        if mcsample!='LcDs':
            bdt = GetBDTcut(BDTcut)
            print('Appling cuts: ', bdt)
            print('... CopyingTree with final preselections ...')
            ot=t.CloneTree(0)
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if t.bdt>BDTcut:
                    if TruthMatch(t):
                        if TruthMatchLambdab(t) and TruthMatchLambdac(t):
                            if TruthMatchMuonPdg(t) and TruthMatchMuonLb(t):
                                ot.Fill()
            of.Write()
            of.Close()
            f.Close()

        else:
            ot = t.CloneTree(0)
            w_2charm = np.zeros(1,dtype=float)

            ot.Branch("w_2charm",w_2charm,"w_2charm/D")
            h_mLcD = r.TH1F('h_mLcD',';m_{L_c D} (MeV^2);',100,3000,6000)
            h_mLcD.SetDirectory(0)

            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if t.bdt>BDTcut:
                    if TruthMatch(t):
                        if TruthMatchLambdab(t) and TruthMatchLambdac(t):
                            if TruthMatchMuonPdg(t) and TruthMatchMuonCharm(t):
                                if SelectMBody(t):
                                    Lc_LV = r.TLorentzVector(t.Lb_TrueHadron_D0_PX,t.Lb_TrueHadron_D0_PY,t.Lb_TrueHadron_D0_PZ, t.Lb_TrueHadron_D0_PE)
                                    D_LV = r.TLorentzVector(t.Lb_TrueHadron_D1_PX,t.Lb_TrueHadron_D1_PY,t.Lb_TrueHadron_D1_PZ, t.Lb_TrueHadron_D1_PE)

                                    mLcD2 = (Lc_LV+D_LV).M2()
                                    h_mLcD.Fill(r.TMath.Sqrt(mLcD2))
                                    mLc = Lc_LV.M()
                                    mD = D_LV.M()
                                    mLb = 5619.6
                                    mK = 493.677
                                    w_2charm[0] = (mLcD2 - rm.Power(mLc+mD,2))/(rm.Power(mLb-mK,2)-rm.Power(mLc+mD,2))
                                    if w_2charm[0]>0:
                                        w_2charm[0] = rm.Sqrt(w_2charm)-1./2
                                    else:
                                        print(' Negative 2charm weight')
                                        w_2charm[0] = -1000.
                                    ot.Fill()
                                if Select2Body(t):
                                    w_2charm[0]=1
                                    ot.Fill()
            of.Write()
            of.Close()
            f.Close()

            c = r.TCanvas('c','mLcD',500,500)
            h_mLcD.Draw()
            c.SaveAs('plots/m_LbD.png')



        print('<<<<< Created file ', ofname)
    return ofname

                
