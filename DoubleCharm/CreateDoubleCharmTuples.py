import ROOT as r

files = {'MagUp':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root','MagDown':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root'}
polarities=['MagUp','MagDown']

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

def Select2Body(tree):
    if abs(int(getattr(tree,'Lb_TrueHadron_D0_ID')))==4122 and int(getattr(tree,'Lb_TrueHadron_D1_ID'))!=0 and int(getattr(tree,'Lb_TrueHadron_D2_ID'))==0:
        flag = True
    else:
        flag = False
    return flag

def SelectMBody(tree):
    if abs(getattr(tree,'Lb_TrueHadron_D0_ID'))==4122 and getattr(tree,'Lb_TrueHadron_D1_ID')!=0 and getattr(tree,'Lb_TrueHadron_D2_ID')!=0:
        flag = True
    else:
        flag = False
    return flag

def TruthMatchMuon(tree):
    if abs(int(getattr(tree,'mu_TRUEID')))==13:
        if abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==431 or abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==421 or  abs(int(getattr(tree,'mu_MC_MOTHER_ID')))==411:
            flag = True
        else:
            flag=False
    else:
        flag = False
    return flag

for polarity in polarities:
    f = r.TFile(files[polarity],'READ')
    t = f.Get('DecayTree')
    of = r.TFile('$FILEDIR/MC/Lb_LcDs_2body_'+polarity+'.root','RECREATE')
    ot = r.TTree('DecayTree','DecayTree')
    ot = t.CloneTree(0)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if TruthMatch(t):
            if TruthMatchLambdab(t) and TruthMatchLambdac(t):
                if Select2Body(t) and TruthMatchMuon(t):
                    ot.Fill()
    of.Write()
    of.Close()
    f.Close()

for polarity in polarities:
    f = r.TFile(files[polarity])
    t = f.Get('DecayTree')
    of = r.TFile('$FILEDIR/MC/Lb_LcDs_mbody_'+polarity+'.root','RECREATE')
    ot = r.TTree('DecayTree','DecayTree')
    ot = t.CloneTree(0)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        if TruthMatch(t):
            if TruthMatchLambdab(t) and TruthMatchLambdac(t):
                if SelectMBody(t) and TruthMatchMuon(t):
                    ot.Fill()
    of.Write()
    of.Close()
    f.Close()




