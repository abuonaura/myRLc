import ROOT as r

def GetTriggerCut(dtype):
    L0 = r.TCut("Lb_L0Global_TIS==1 || Lc_L0HadronDecision_TOS == 1")
    HLT1 = r.TCut("Lc_Hlt1TrackMVADecision_TOS==1 || Lc_Hlt1TwoTrackMVADecision_TOS==1")
    if dtype =='Data' or dtype == 'DataSS' or dtype == 'MC':
        HLT2 = r.TCut("Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==1")
    else:
        HLT2 = r.TCut("Lb_Hlt2XcMuXForTauB2XcFakeMuDecision_TOS==1")
    Trigger = r.TCut("("+L0.GetTitle() + ') && (' + HLT1.GetTitle() + ') && (' +  HLT2.GetTitle()+")") 
    return Trigger


def GetStrippingPIDcut():
    pidcut = r.TCut("pi_PIDK_corr<2 && K_PIDK_corr>4. && p_PIDp_corr>0.")
    return pidcut

def ReduceTree(fname, tname, nfname, ntname, dtype):
    f = r.TFile(fname, 'READ')
    if not f:
        print('File not found!')
    t = f.Get(tname)
    if not t:
        print('Tree not found!')
    if dtype!='MC':
        tl = f.Get('GetIntegratedLuminosity/LumiTuple')

    of = r.TFile(nfname, 'RECREATE')
    ot = r.TTree(ntname, ntname)
    if dtype!='MC':
        print('... Adding Luminosity Tree ...')
        otl = r.TTree('LumiTuple','LumiTuple')
        otl = tl.CloneTree()
        
    
    cTrigger = GetTriggerCut(dtype)
    cMass = r.TCut("Lc_M>2230 && Lc_M<2330")
    if dtype!='MC':
        cStep1 = r.TCut("("+cTrigger.GetTitle()+ ') && (' + cMass.GetTitle()+")")
    else:
        pidcut = GetStrippingPIDcut()
        cStep1 = r.TCut("("+cTrigger.GetTitle()+ ') && (' + cMass.GetTitle()+") && ("+ pidcut.GetTitle()+")")

    print ('Applied cuts: ', cStep1)

    ot = t.CopyTree(cStep1.GetTitle())
    of.Write()
    of.Close()
    f.Close()
    return cTrigger

def TruthMatch(tree):
    if int(getattr(tree,'Lc_BKGCAT'))<30 and int(getattr(tree,'Lb_BKGCAT'))<50:
        flag = True
    else:
        flag = False
    return flag

def TruthMatchLambdab(tree):
    if abs(int(getattr(tree,'Lb_TRUEID')))==5122:
        flag=True
    else:
        flag=False
    return flag

def TruthMatchLambdac(tree):
    if abs(int(getattr(tree,'Lc_TRUEID')))==4122:
        flag=True
    else:
        flag=False
    return flag

def TruthMatchMuonPdg(tree):
    if abs(int(getattr(tree,'mu_TRUEID')))==13:
        flag=True
    else:
        flag = False
    return flag

def TruthMatchMuonCharm(tree):
    if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==431 or abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==421 or  abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==411:
        flag = True
    else:
        flag=False
    return flag

def TruthMatchMuonLb(tree):
    if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==5122:
        flag = True
    else:
        flag=False
    return flag



