import ROOT as r

files = {'MagUp':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root','MagDown':'$FILEDIR/MC/Lb_LcDs_MagUp_PID_reduced_preselected.root'}
polarities=['MagUp','MagDown']

def GetTotalNumberMCevts():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.GetEntries()
    print('Total number of MC events before applying ANY selection: ', ntotMC)
    return ntotMC

def GetTotalNumber2body():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  t.Lb_TrueHadron_D2_ID==0:
                ntot_2body+=1
    print('Total number of 2body MC events before applying ANY selection: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  t.Lb_TrueHadron_D2_ID!=0:
                ntot_mbody+=1
    print('Total number of mbody MC events before applying ANY selection: ', ntot_mbody)
    return ntot_mbody

def TruthMatch(tree):
    if getattr(tree,'Lc_BKGCAT')<30 and getattr(tree,'Lb_BKGCAT')<50:
        flag = True
    else:
        flag = False
    return flag

def GetTotalNumberMCevts_TMatch():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if TruthMatch(t)==True:
                ntotMC+=1
    print('Total number of MC events before applying ANY selection: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if TruthMatch(t)==True and t.Lb_TrueHadron_D2_ID==0:
                ntot_2body+=1
    print('Total number of 2body MC events before applying ANY selection: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  TruthMatch(t)==True and t.Lb_TrueHadron_D2_ID!=0:
                ntot_mbody+=1
    print('Total number of mbody MC events before applying ANY selection: ', ntot_mbody)
    return ntot_mbody




GetTotalNumberMCevts()
GetTotalNumber2body()
GetTotalNumberMbody()
print()
GetTotalNumberMCevts_TMatch()
GetTotalNumber2body_TMatch()
GetTotalNumberMbody_TMatch()
