import ROOT as r
from ROOT import TMath as tm
from root_pandas import read_root
from array import array
import sys, os
from datetime import datetime
from Trigger import *
from TruthMatch import *
from PIDcuts import *
from SplitSamples import *
sys.path.append('../PIDGen_PIDCalib_MVA/')
from Add_MVA import AddBDTinfo

BDTcut=0.7
filedir = '/disk/lhcb_data2/RLcMuonic2016/'

def LcMussCut(t,evt):
    t.GetEntry(evt)
    if t.Lc_M>2230 and t.Lc_M<2330:
        return 1
    else:
        return 0

def GetPIDCalibCuts(t, evt):
    t.GetEntry(evt)
    if (t.nTracks<700 and t.nTracks>0) and (t.mu_P>0 and t.mu_P<200000) and (t.mu_PT>0 and t.mu_PT<60000) and (t.pi_P>0 and t.pi_P<200000) and (t.pi_PT>0 and t.pi_PT<60000) and (t.p_P>0 and t.p_P<200000) and (t.p_PT>0 and t.p_PT<60000) and (t.K_P>0 and t.K_P<200000) and (t.K_PT>0 and t.K_PT<60000):
        return 1
    else:
        return 0

def GetFinalPreselectionData(trigger,LcMass,PIDcut,PIDCalibCut):
    if trigger==1 and LcMass==1 and PIDcut==1 and PIDCalibCut==1:
        return 1
    else:
        return 0


def CheckSurvivalBDT(bdt,BDTcut):
    if bdt>BDTcut: 
        return 1
    else:
        return 0

def GetFinalSelectionData(finalPreselection,passBDT,RemoveDDst):
    if finalPreselection==1 and passBDT==1 and RemoveDDst==1:
        return 1
    else: 
        return 0

def GetRunEventNumber(tree,i):
    tree.GetEntry(i)
    return tree.runNumber,tree.eventNumber


def CreateSelectionTreeData(ifname,tree,dtype,polarity,nevts):
    print(ifname)
    bdtfname = ifname[0:-5]+'_MVA.root'
    if os.path.isfile(bdtfname):
        print('BDT file already created')
    else:
        print()
        print('>>>   Creating file with BDT variable')
        print()
        AddBDTinfo(ifname, 'tupleout/DecayTree', bdtfname, 'Data',pickled_model_path = '../PIDGen_PIDCalib_MVA/xgb_reg.pkl')
    df_bdt = read_root(bdtfname,'DecayTree')
    print('HERE')

    ofname = ifname[0:-5]+'_Selections_'+datetime.now().strftime("%d-%m-%Y_%I-%M-%S_%p")+'.root'
    of = r.TFile(ofname,'RECREATE')
    ot = r.TTree('DecayTree','DecayTree')
    L0 = array( 'i', [ 0 ] )
    HLT1 = array( 'i', [ 0 ] )
    HLT2 = array( 'i', [ 0 ] )
    TFinal = array('i',[0])
    LcMass = array('i',[0]) 
    PIDcut = array('i',[0])
    PIDCalibCut = array('i',[0])
    TOTpreselection = array('i',[0])
    PassBDT = array('i',[0])
    RemoveDDst = array('i',[0])
    FinalSelection = array('i',[0])
    isIsolated = array('i',[0])
    isKenriched = array('i',[0])
    ot.Branch('L0',L0,'L0/I')
    ot.Branch('HLT1',HLT1,'HLT1/I')
    ot.Branch('HLT2',HLT2,'HLT2/I')
    ot.Branch('Trigger_Final',TFinal,'Trigger_Final/I')
    ot.Branch('LcMass',LcMass,'LcMass/I')
    ot.Branch('PIDcut',PIDcut,'PIDcut/I')
    ot.Branch('PIDCalibCut',PIDCalibCut,'PIDCalibCut/I')
    ot.Branch('TOTpreselection',TOTpreselection,'TOTpreselection/I')
    ot.Branch('PassBDT',PassBDT,'PassBDT/I')
    ot.Branch('RemoveDDst',RemoveDDst,'RemoveDDst/I')
    ot.Branch('FinalSelection',FinalSelection,'FinalSelection/I')
    ot.Branch('isIsolated',isIsolated,'isIsolated/I')
    ot.Branch('isKenriched',isKenriched,'isKenriched/I')
    

    t.SetBranchStatus('*',0)
    varsON = ['Lb_L0Global_TIS','Lb_L0HadronDecision_TOS','Lc_Hlt1TrackMVADecision_TOS','Lc_Hlt1TwoTrackMVADecision_TOS','Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS','Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS','Lc_M','p_ProbNNp','p_ProbNNk','mu_PID*','*_P','*_PT','nTracks','runNumber','eventNumber','Lb_ISOLATION_*','mu_PX','mu_PY','mu_PZ','mu_ID','Lc_PX','Lc_PY','Lc_PZ']
    for var in varsON:
        t.SetBranchStatus(var,1)

    for i in range(nevts):
        if i%500000==0:
            print(i)
        L0[0] = L0TriggerData(tree,i)
        HLT1[0] = HLT1TriggerData(tree,i)
        HLT2[0] = HLT2TriggerData(tree,i,dtype)
        TFinal[0] = GetTriggerData(L0[0],HLT1[0],HLT2[0])
        LcMass[0]= LcMussCut(tree,i)
        PIDcut[0] = GetPIDCutData(tree,i,dtype)
        PIDCalibCut[0] = GetPIDCalibCuts(tree, i)
        TOTpreselection[0] = GetFinalPreselectionData(TFinal[0],LcMass[0],PIDcut[0],PIDCalibCut[0])
        runN, evtN = GetRunEventNumber(tree,i)
        #bdt = df_bdt.loc[(df_bdt['eventNumber']==evtN)&(df_bdt['runNumber']==runN),'bdt']
        runN_2, evtN_2 = df_bdt.loc[i,'runNumber'], df_bdt.loc[i,'eventNumber']
        if runN!=runN_2 or evtN!=evtN_2:
            print('ERROR!!!!')
        #print(runN, runN_2, evtN, evtN_2)
        bdt = df_bdt.loc[i,'bdt']
        #print(bdt)
        PassBDT[0] = CheckSurvivalBDT(bdt,BDTcut)
        RemoveDDst[0] = CheckPProbNNcut(tree,i)
        FinalSelection[0] = GetFinalSelectionData(TOTpreselection[0],PassBDT[0],RemoveDDst[0]) 
        isKenriched[0] = CheckIfIsKenriched(tree,i)
        isIsolated[0] = CheckIfIsIsolated(tree,i)
        ot.Fill()

    of.Write()
    of.Close()


if __name__ == "__main__":
    dtype = ['Data','DataSS','FakeMu','FakeMuSS']
    #dtype = ['DataSS']
    polarities=['MagUp','MagDown']
    #polarities=['MagUp']
    for polarity in polarities:
        for dt in dtype:
            ifname = filedir+'Data/Lb_'+dt+'_'+polarity+'.root'
            f = r.TFile(ifname,'READ')
            t = f.Get('tupleout/DecayTree')
            nevts = t.GetEntries()
            #nevts = 1000000
            CreateSelectionTreeData(ifname,t,dt,polarity,nevts)

