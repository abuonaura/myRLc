import ROOT as r

#folder = '/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/'
folder = '/Users/annarita/cernbox/LHCbDatasets/MC_full_TrueIsoInfo/'
files = {'MagUp':folder + 'Lb_LcDs_MagUp_full.root','MagDown':folder + 'Lb_LcDs_MagDown_full.root'}
preselfiles = {'MagUp': folder + 'Lb_LcDs_MagUp_full_preselectionVars.root', 'MagDown': folder + 'Lb_LcDs_MagDown_full_preselectionVars.root'}
polarities=['MagUp','MagDown']



def GetTotalNumberMCevts():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntotMC += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lc_M>0)','goff')
    print('Total number of MC events before applying ANY selection (weighted): ', ntotMC)
    return ntotMC

def GetTotalNumber2body():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntot_2body += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID==0)','goff')
    print('Total number of 2body MC events before applying ANY selection (weighted): ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntot_mbody += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID!=0)','goff')
    print('Total number of mbody MC events before applying ANY selection (weighted): ', ntot_mbody)
    return ntot_mbody

def GetTotalNumberMCevts_TMatch():
    ntotMC=0 #total number  of MC events
    ntot_2body = 0
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntotMC += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(TruthMatch==1)','goff')
        ntot_2body+=t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID==0&&TruthMatch==1)','goff')
        ntot_mbody += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID!=0&&TruthMatch==1)','goff')
    print('Total number of MC events after Truth Matching: ', ntotMC)
    print('Total number of 2body MC evts after Truth Matching: ', ntot_2body)
    print('Total number of mbody MC evts after Truth Matching: ', ntot_mbody)
    return

def GetTotalNumberMCevts_Trigger():
    ntotMC=0 #total number  of MC events
    ntot_2body = 0
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntotMC += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(TruthMatch==1&&Trigger==1)','goff')
        ntot_2body+=t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID==0&&TruthMatch==1&&Trigger==1)','goff')
        ntot_mbody += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID!=0&&TruthMatch==1&&Trigger==1)','goff')
    print('Total number of MC events after Truth Matching + Trigger: ', ntotMC)
    print('Total number of 2body MC evts after Truth Matching + Trigger: ', ntot_2body)
    print('Total number of mbody MC evts after Truth Matching + Trigger: ', ntot_mbody)
    return

def GetTotalNumberMCevts_FullPreselection():
    ntotMC=0 #total number  of MC events
    ntot_2body = 0
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        ntotMC += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(FinalSel==1)','goff')
        ntot_2body+=t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID==0&&FinalSel==1)','goff')
        ntot_mbody += t.Draw('Lc_M','Event_PIDCalibEffWeight*w_LbCorr*(Lb_TrueHadron_D0_ID==4122&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID!=0&&FinalSel==1)','goff')
    print('Total number of MC events after full preselection: ', ntotMC)
    print('Total number of 2body MC evts after full preselection: ', ntot_2body)
    print('Total number of mbody MC evts after full preselection: ', ntot_mbody)
    return


def GetTotalNumberMCevts_TMatch_Iso():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of MC events after applying Truth Matching + Isolation: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch_Iso():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID==0&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of 2body MC events after applying Truth Matching + Isolation: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch_Iso():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID!=0&&Lb_ISOLATION_BDT<0.35','goff')
    print('Total number of mbody MC events after applying Truth Matching + Isolation: ', ntot_mbody)
    return ntot_mbody

def GetTotalNumberMCevts_TMatch_Kenriched():
    ntotMC=0 #total number  of MC events
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntotMC += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50','goff')
    print('Total number of MC events after applying Truth Matching +  Kenriched selection: ', ntotMC)
    return ntotMC

def GetTotalNumber2body_TMatch_Kenriched():
    ntot_2body = 0
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_2body += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID==0','goff')
    print('Total number of 2body MC events after applying Truth Matching +  Kenriched selection: ', ntot_2body)
    return ntot_2body

def GetTotalNumberMbody_TMatch_Kenriched():
    ntot_mbody = 0
    for polarity in polarities:
        f = r.TFile(csamples[polarity],'READ')
        t = f.Get('DecayTree')
        ntot_mbody += t.Draw('Lc_M','Lc_BKGCAT<30 && Lb_BKGCAT<50&&Lb_TrueHadron_D2_ID!=0','goff')
    print('Total number of mbody MC events after applying Truth Matching + Kenriched selection : ', ntot_mbody)
    return ntot_mbody

def FillHistoD1(h,pdg):
    if r.TMath.Abs(pdg)==411:
        h.Fill(1)
    if r.TMath.Abs(pdg)==421:
        h.Fill(2)
    if r.TMath.Abs(pdg)==413:
        h.Fill(3)
    if r.TMath.Abs(pdg)==423:
        h.Fill(4)
    if r.TMath.Abs(pdg)==415:
        h.Fill(5)
    if r.TMath.Abs(pdg)==425:
        h.Fill(6)
    if r.TMath.Abs(pdg)==431:
        h.Fill(7)
    if r.TMath.Abs(pdg)==10431:
        h.Fill(8)
    if r.TMath.Abs(pdg)==433:
        h.Fill(9)
    if r.TMath.Abs(pdg)==435:
        h.Fill(10)
    h.SetDirectory(0)
    return h

def FillHistoD2(h,pdg):
    if r.TMath.Abs(pdg)==311:
        h.Fill(1)
    if r.TMath.Abs(pdg)==321:
        h.Fill(2)
    if r.TMath.Abs(pdg)==10311:
        h.Fill(3)
    h.SetDirectory(0)
    return h

def PutTogetherPolarityHistos(h):
    h_new = h['MagUp']
    h_new.Add(h['MagDown'])
    h_new.SetDirectory(0)
    return h_new

def DisplayLbTrueHadronD0ID(cut=''):
    hD0ID = {polarity: r.TH1F('h_trueD0ID_'+polarity,'',100,0,5000) for polarity in polarities}
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        name = polarity
        if cut!='':
            t.Draw('Lb_TrueHadron_D0_ID>>D0ID_'+name+'(100,0,5000)','Event_PIDCalibEffWeight*w_LbCorr*('+cut+')')
        else:
            t.Draw('Lb_TrueHadron_D0_ID>>D0ID_'+name+'(100,0,5000)','Event_PIDCalibEffWeight*w_LbCorr*(Lc_M>0)')

        hD0ID[polarity] = r.gPad.GetPrimitive('D0ID_'+name)
        hD0ID[polarity].SetDirectory(0)

    hD0ID_1 = PutTogetherPolarityHistos(hD0ID)
    return hD0ID_1
        
def DisplayLbTrueHadronD1ID(cut=''):
    hD1ID = {polarity: r.TH1F('h_trueD1ID_'+polarity,'',500,0,500) for polarity in polarities}
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        name = polarity
        if cut!='':
            t.Draw('Lb_TrueHadron_D1_ID>>D1ID_'+name+'(500,0,500)','Event_PIDCalibEffWeight*w_LbCorr*('+cut+')')
        else:
            t.Draw('Lb_TrueHadron_D1_ID>>D1ID_'+name+'(500,0,500)','Event_PIDCalibEffWeight*w_LbCorr*(Lc_M>0)')

        hD1ID[polarity] = r.gPad.GetPrimitive('D1ID_'+name)
        hD1ID[polarity].SetDirectory(0)

    hD1ID_1 = PutTogetherPolarityHistos(hD1ID)
    return hD1ID_1

def DisplayLbTrueHadronD2ID(cut=''):
    hD2ID = {polarity: r.TH1F('h_trueD2ID_'+polarity,'',500,0,500) for polarity in polarities}
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        name = polarity
        if cut!='':
            t.Draw('Lb_TrueHadron_D2_ID>>D2ID_'+name+'(500,0,500)','Event_PIDCalibEffWeight*w_LbCorr*('+cut+')')
        else:
            t.Draw('Lb_TrueHadron_D2_ID>>D2ID_'+name+'(500,0,500)','Event_PIDCalibEffWeight*w_LbCorr*(Lc_M>0)')

        hD2ID[polarity] = r.gPad.GetPrimitive('D2ID_'+name)
        hD2ID[polarity].SetDirectory(0)

    hD2ID_1 = PutTogetherPolarityHistos(hD2ID)
    return hD2ID_1


'''

c = r.TCanvas('c','',1500,500)
c.Divide(3,1)
c.cd(1)
h0 =  DisplayLbTrueHadronD0ID('')
h0.Draw('hist')
c.cd(2)
h0a =  DisplayLbTrueHadronD0ID('TruthMatch==1&&Trigger==1')
h0a.Draw('hist')
c.cd(3)
h0b =  DisplayLbTrueHadronD0ID('FinalSel==1')
h0b.Draw('hist')

c1 = r.TCanvas('c1','',1500,500)
c1.Divide(3,1)
c1.cd(1)
h1 =  DisplayLbTrueHadronD1ID('')
h1.Draw('hist')
c1.cd(2)
h1a =  DisplayLbTrueHadronD1ID('TruthMatch==1&&Trigger==1')
h1a.Draw('hist')
c1.cd(3)
h1b =  DisplayLbTrueHadronD1ID('FinalSel==1')
h1b.Draw('hist')


c2 = r.TCanvas('c2','',1500,500)
c2.Divide(3,1)
c2.cd(1)
h2 =  DisplayLbTrueHadronD2ID('')
h2.Draw('hist')
c2.cd(2)
h2a =  DisplayLbTrueHadronD2ID('TruthMatch==1&&Trigger==1')
h2a.Draw('hist')
c2.cd(3)
h2b =  DisplayLbTrueHadronD2ID('FinalSel==1')
h2b.Draw('hist')
'''

def PrintDecays(nentries):
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
            nentries = t.GetEntries()
        for i in range(nentries):
            t.GetEntry(i)
            if t.FinalSel!=1:
                continue
            else:
                print('-----------------------------------------------')
                print('Evt: '+ str(i))
                print('LbID: '+str(t.Lb_TRUEID)+' D0: '+str(t.Lb_TrueHadron_D0_ID)+' D1: '+str(t.Lb_TrueHadron_D1_ID)+ ' D2: '+str(t.Lb_TrueHadron_D2_ID))
                print('- '+str(t.Lb_TrueHadron_D0_ID)+' -> '+str(t.Lb_TrueHadron_D0_GD0_ID)+' '+str(t.Lb_TrueHadron_D0_GD1_ID)+ ' '+str(t.Lb_TrueHadron_D0_GD2_ID))
                print('- '+str(t.Lb_TrueHadron_D1_ID)+' -> '+str(t.Lb_TrueHadron_D1_GD0_ID)+' '+str(t.Lb_TrueHadron_D1_GD1_ID)+ ' '+str(t.Lb_TrueHadron_D1_GD2_ID))
                if(t.Lb_TrueHadron_D2_ID!=0):
                    print('- '+str(t.Lb_TrueHadron_D2_ID)+' -> '+str(t.Lb_TrueHadron_D2_GD0_ID)+' '+str(t.Lb_TrueHadron_D2_GD1_ID)+ ' '+str(t.Lb_TrueHadron_D2_GD2_ID))
                PrintLbIsoParticles(t, i)
    return

def PrintLbIsoParticles(t, i):
    t.GetEntry(i)
    print(' 1. BDT:' +str(t.Lb_ISOLATION_BDT)+' TruePID: '+str(t.Lb_ISOLATION_TruePID)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID))
    print(' 2. BDT:' +str(t.Lb_ISOLATION_BDT2)+' TruePID: '+str(t.Lb_ISOLATION_TruePID2)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID2)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID2))
    print(' 3. BDT:' +str(t.Lb_ISOLATION_BDT3)+' TruePID: '+str(t.Lb_ISOLATION_TruePID3)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID3)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID3))
    return

KIDs = [321,311,313,310,323]

def CountTrueKaonsInDecays(nentries):
    nK1ryDecay,nK1stDaughterDecay,nK2ndDaughterDecay = 0,0,0
    for polarity in polarities:
        print(polarity)
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
            nentries = t.GetEntries()
            print(nentries)
        for i in range(nentries):
            t.GetEntry(i)
            if i%100000==0:
                print(i)
            if t.FinalSel==1:
                if (abs(t.Lb_TrueHadron_D1_ID) in KIDs) or (abs(t.Lb_TrueHadron_D2_ID) in KIDs):
                    nK1ryDecay+=1
                if (abs(t.Lb_TrueHadron_D1_GD0_ID) in KIDs) or (abs(t.Lb_TrueHadron_D1_GD1_ID) in KIDs) or (abs(t.Lb_TrueHadron_D1_GD2_ID) in KIDs):
                    nK1stDaughterDecay+=1
                if (abs(t.Lb_TrueHadron_D2_GD0_ID) in KIDs) or (abs(t.Lb_TrueHadron_D2_GD1_ID) in KIDs) or (abs(t.Lb_TrueHadron_D2_GD2_ID) in KIDs):
                    nK2ndDaughterDecay+=1
    
    print('Number of evts with true K from Lb decay: ', nK1ryDecay)
    print('Number of evts with true K from decay of the 1st Lb daughter: ',nK1stDaughterDecay)
    print('Number of evts with true K from decay of the 2nd Lb daughter: ',nK2ndDaughterDecay)
    return

#PrintDecays(100)
CountTrueKaonsInDecays(-1)
#PrintLbIsoDecays(1000)

def PrintDecays(t,entry):
    t.GetEntry(entry)
    print('LbID: '+str(t.Lb_TRUEID)+' D0: '+str(t.Lb_TrueHadron_D0_ID)+' D1: '+str(t.Lb_TrueHadron_D1_ID)+ ' D2: '+str(t.Lb_TrueHadron_D2_ID))
    print(' 1. BDT:' +str(t.Lb_ISOLATION_BDT)+' TruePID: '+str(t.Lb_ISOLATION_TruePID)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID)) 
    print(' 2. BDT:' +str(t.Lb_ISOLATION_BDT2)+' TruePID: '+str(t.Lb_ISOLATION_TruePID2)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID2)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID2)) 
    print(' 3. BDT:' +str(t.Lb_ISOLATION_BDT3)+' TruePID: '+str(t.Lb_ISOLATION_TruePID3)+ ' _Mother: '+str(t.Lb_ISOLATION_TrueMotherPID3)+' _Grandma: '+str(t.Lb_ISOLATION_TrueGrandmotherPID3)) 
    return


def AnalysisIsoParticles(nentries):
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
            nentries = t.GetEntries()
        for i in range(nentries):
            t.GetEntry(i)
            if t.FinalSel!=1:
                continue
            else:
                #Lambda_b daughters
                LbDaughters = [t.Lb_TrueHadron_D0_ID, t.Lb_TrueHadron_D1_ID, t.Lb_TrueHadron_D2_ID]
                #particles from isolation algorithm
                isoparticles = [t.Lb_ISOLATION_TruePID, t.Lb_ISOLATION_TruePID2, t.Lb_ISOLATION_TruePID3]
                #mothers of particles from isolation algorithm
                misoparticles = [t.Lb_ISOLATION_TrueMotherPID, t.Lb_ISOLATION_TrueMotherPID2, t.Lb_ISOLATION_TrueMotherPID3]
                #grandmothers of particles from isolation algorithm
                gmisoparticles = [t.Lb_ISOLATION_TrueGrandmotherPID, t.Lb_ISOLATION_TrueGrandmotherPID2, t.Lb_ISOLATION_TrueGrandmotherPID3]
                evtfound_m, evtfound_gm = 0,0
                print('---------------------------------------')
                for j,isop in enumerate(isoparticles):
                    if isop!=0:
                        if misoparticles[j] in LbDaughters and misoparticles[j]!=0:
                            print(isop, misoparticles[j])
                            if evtfound_m==0:
                                PrintDecays(t,i)
                            evtfound_m=1
                        if evtfound_m==0:
                            if gmisoparticles[j] in LbDaughters and gmisoparticles[j]!=0:
                                print(isop, misoparticles[j], gmisoparticles[j])
                                if evtfound_gm==0:
                                    PrintDecays(t,i)
                                evtfound_gm=1
                print('---------------------------------------')
    return

#AnalysisIsoParticles(1000)

def CountIsoParicleisLbDaughter(nentries,ndaughter):
    count, count_1, count_2, count_3 = 0, 0, 0, 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
             nentries = t.GetEntries()
        for i in range(nentries):
            t.GetEntry(i)
            if t.FinalSel!=1:
                continue
            else:
                #Lambda_b daughters
                LbDaughters = [t.Lb_TrueHadron_D0_ID, t.Lb_TrueHadron_D1_ID, t.Lb_TrueHadron_D2_ID]
                #particles from isolation algorithm
                isoparticles = [t.Lb_ISOLATION_TruePID, t.Lb_ISOLATION_TruePID2, t.Lb_ISOLATION_TruePID3]
                #mothers of particles from isolation algorithm
                misoparticles = [t.Lb_ISOLATION_TrueMotherPID, t.Lb_ISOLATION_TrueMotherPID2, t.Lb_ISOLATION_TrueMotherPID3]
                for j,isop in enumerate(isoparticles):
                    if isop!=0:
                        if isop==LbDaughters[ndaughter] and misoparticles[j]==t.Lb_TRUEID:
                            count+=1
                            if j==0:
                                count_1+=1
                            if j==1:
                                count_2+=1
                            if j==2:
                                count_3+=1
        print('Polarity :',polarity)
        print('Number of evts where one isolated particle is Lb daughter '+str(ndaughter)+' '+str(count))
        print('Number of evts where 1st isolated particle is Lb daughter '+str(ndaughter)+' '+str(count_1))
        print('Number of evts where 2nd isolated particle is Lb daughter '+str(ndaughter)+' '+str(count_2))
        print('Number of evts where 3rd isolated particle is Lb daughter '+str(ndaughter)+' '+str(count_3))
    return


def CountIsoParicleisLbGdaughter(nentries,ndaughter):
    count, count_1, count_2, count_3 = 0, 0, 0, 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
             nentries = t.GetEntries()
        for i in range(nentries):
            t.GetEntry(i)
            if t.FinalSel!=1:
                continue
            else:
                #Lambda_b daughters
                LbDaughters = [t.Lb_TrueHadron_D0_ID, t.Lb_TrueHadron_D1_ID, t.Lb_TrueHadron_D2_ID]
                #particles from isolation algorithm
                isoparticles = [t.Lb_ISOLATION_TruePID, t.Lb_ISOLATION_TruePID2, t.Lb_ISOLATION_TruePID3]
                #mothers of particles from isolation algorithm
                misoparticles = [t.Lb_ISOLATION_TrueMotherPID, t.Lb_ISOLATION_TrueMotherPID2, t.Lb_ISOLATION_TrueMotherPID3]
                #grandmothers of particles from isolation algorithm
                gmisoparticles = [t.Lb_ISOLATION_TrueGrandmotherPID, t.Lb_ISOLATION_TrueGrandmotherPID2, t.Lb_ISOLATION_TrueGrandmotherPID3]
                for j,isop in enumerate(isoparticles):
                    if isop!=0:
                        if misoparticles[j]==LbDaughters[ndaughter] and gmisoparticles[j]==t.Lb_TRUEID:
                            count+=1
                            if j==0:
                                count_1+=1
                            if j==1:
                                count_2+=1
                            if j==2:
                                count_3+=1
        print('Polarity :',polarity)
        print('Number of evts where one isolated particle is daughter of Lb Daughter '+str(ndaughter)+' '+str(count))
        print('Number of evts where 1st isolated particle is daughter of Lb Daughter '+str(ndaughter)+' '+str(count_1))
        print('Number of evts where 2nd isolated particle is daughter of Lb Daughter '+str(ndaughter)+' '+str(count_2))
        print('Number of evts where 3rd isolated particle is daughter of Lb Daughter '+str(ndaughter)+' '+str(count_3))
    return

def CountIsoParicleisLbGGdaughter(nentries,ndaughter):
    count, count_1, count_2, count_3 = 0, 0, 0, 0
    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('tupleout/DecayTree')
        fpresel = r.TFile(preselfiles[polarity],'READ')
        tpresel = fpresel.Get('DecayTree')
        t.AddFriend(tpresel)
        if nentries==-1:
             nentries = t.GetEntries()
        for i in range(nentries):
            t.GetEntry(i)
            if t.FinalSel!=1:
                continue
            else:
                #Lambda_b daughters
                LbDaughters = [t.Lb_TrueHadron_D0_ID, t.Lb_TrueHadron_D1_ID, t.Lb_TrueHadron_D2_ID]
                #particles from isolation algorithm
                isoparticles = [t.Lb_ISOLATION_TruePID, t.Lb_ISOLATION_TruePID2, t.Lb_ISOLATION_TruePID3]
                #mothers of particles from isolation algorithm
                misoparticles = [t.Lb_ISOLATION_TrueMotherPID, t.Lb_ISOLATION_TrueMotherPID2, t.Lb_ISOLATION_TrueMotherPID3]
                #grandmothers of particles from isolation algorithm
                gmisoparticles = [t.Lb_ISOLATION_TrueGrandmotherPID, t.Lb_ISOLATION_TrueGrandmotherPID2, t.Lb_ISOLATION_TrueGrandmotherPID3]
                for j,isop in enumerate(isoparticles):
                    if isop!=0:
                        if misoparticles[j]!=0 and gmisoparticles[j]==LbDaughters[ndaughter]:
                            count+=1
                            if j==0:
                                count_1+=1
                            if j==1:
                                count_2+=1
                            if j==2:
                                count_3+=1
        print('Polarity :',polarity)
        print('Number of evts where one isolated particle is (maybe) granddaughter of Lb Daughter '+str(ndaughter)+' '+str(count))
        print('Number of evts where 1st isolated particle is (maybe) granddaughter of Lb Daughter '+str(ndaughter)+' '+str(count_1))
        print('Number of evts where 2nd isolated particle is (maybe) granddaughter of Lb Daughter '+str(ndaughter)+' '+str(count_2))
        print('Number of evts where 3rd isolated particle is (maybe) granddaughter of Lb Daughter '+str(ndaughter)+' '+str(count_3))
    return

'''
CountIsoParicleisLbDaughter(10000,0)
CountIsoParicleisLbGdaughter(10000,0)
CountIsoParicleisLbGGdaughter(10000,0)
print()

CountIsoParicleisLbDaughter(10000,1)
CountIsoParicleisLbGdaughter(10000,1)
CountIsoParicleisLbGGdaughter(10000,1)
print()
CountIsoParicleisLbDaughter(10000,2)
CountIsoParicleisLbGdaughter(10000,2)
CountIsoParicleisLbGGdaughter(10000,2)
'''

def Check_mbody():
    h_BDT = r.TH1F('h_BDT','Lb_ISOLATION_BDT',50,-2,1)
    h_BDT_sel = r.TH1F('h_BDT_sel','Lb_ISOLATION_BDT',50,-2,1)
    h_BDT2 = r.TH1F('h_BDT2','Lb_ISOLATION_BDT2',50,-2,1)
    h_BDT2_sel = r.TH1F('h_BDT2_sel','Lb_ISOLATION_BDT2',50,-2,1)
    h_TrueD1 = r.TH1F('h_TrueD1','Lb_TrueHadron_D1_ID',10,1,11)
    h_TrueD1_1= r.TH1F('h_TrueD1_1','Lb_TrueHadron_D1_ID',10,1,11) #PASSING ISO_BDT>0.35 AND ISO_BDT2>0.2
    D1particles = ['D^{+}','D^{0}','D^{*}(2010)^{+}','D^{*}(2007)^{0}','D^{*}_{2}(2460)^{+}','D^{*}_{2}(2460)^{0}','D^{+}_{s}','D^{*}_{s0}(2317)^{+}','D^{*+}_{s}','D^{*}_{s}(2573)^{+}']
    for n, d1p in enumerate(D1particles):
        print(n, d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1_1.GetXaxis().SetBinLabel(n+1,d1p)
        h_TrueD1_1.GetXaxis().LabelsOption('v')
        h_TrueD1_1.GetXaxis().LabelsOption('v')
    h_TrueD2 = r.TH1F('h_TrueD2','Lb_TrueHadron_D2_ID',3,1,4)
    h_TrueD2_1 = r.TH1F('h_TrueD2_1','Lb_TrueHadron_D2_ID',3,1,4)  #PASSING ISO_BDT>0.35 AND ISO_BDT2>0.2
    D2particles = ['K^{0}','K','K^{0*}']
    for n, d2p in enumerate(D2particles):
        print(n, d2p)
        h_TrueD2.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2.GetXaxis().LabelsOption('v')
        h_TrueD2_1.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2_1.GetXaxis().SetBinLabel(n+1,d2p)
        h_TrueD2_1.GetXaxis().LabelsOption('v')
    h_PIDK = r.TH1F('h_PIDK','Lb_ISOLATION_PIDK',50,-200,200)
    h_PIDK_1 = r.TH1F('h_PIDK_1','Lb_ISOLATION_PIDK',50,-200,200) #PASSING BOTH ISO BDT CUTS 
    h_PIDK_2 = r.TH1F('h_PIDK_2','Lb_ISOLATION_PIDK',50,-200,200) #PASSING BOTH ISO BDT CUTS + PIDK>4 (= K)
    h_PIDK2 = r.TH1F('h_PIDK2','Lb_ISOLATION_PIDK2',50,-200,200)
    h_PIDK2_1 = r.TH1F('h_PIDK2_1','Lb_ISOLATION_PIDK2',50,-200,200) #PASSING BOTH ISO BDT CUTS 
    h_PIDK2_2 = r.TH1F('h_PIDK2_2','Lb_ISOLATION_PIDK2',50,-200,200) #PASSING BOTH ISO BDT CUTS + PIDK>4 (= K)

    for polarity in polarities:
        f = r.TFile(files[polarity],'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if  t.Lb_TrueHadron_D2_ID!=0 and t.Lc_BKGCAT<30 and t.Lb_BKGCAT<50:
                h_BDT.Fill(t.Lb_ISOLATION_BDT)
                h_BDT2.Fill(t.Lb_ISOLATION_BDT2)
                h_TrueD1=FillHistoD1(h_TrueD1,t.Lb_TrueHadron_D1_ID)
                h_TrueD2=FillHistoD2(h_TrueD2,t.Lb_TrueHadron_D2_ID)
                h_PIDK.Fill(t.Lb_ISOLATION_PIDK)
                h_PIDK2.Fill(t.Lb_ISOLATION_PIDK2)
                if t.Lb_ISOLATION_BDT>0.35:
                    h_BDT_sel.Fill(t.Lb_ISOLATION_BDT)
                if t.Lb_ISOLATION_BDT2>0.2:
                    h_BDT2_sel.Fill(t.Lb_ISOLATION_BDT2)
                if t.Lb_ISOLATION_BDT>0.35 and t.Lb_ISOLATION_BDT2>0.2:
                    h_TrueD1_1=FillHistoD1(h_TrueD1_1,t.Lb_TrueHadron_D1_ID)
                    h_TrueD2_1=FillHistoD2(h_TrueD2_1,t.Lb_TrueHadron_D2_ID)
                    h_PIDK_1.Fill(t.Lb_ISOLATION_PIDK)
                    h_PIDK2_1.Fill(t.Lb_ISOLATION_PIDK2)
                    if t.Lb_ISOLATION_PIDK>4:
                        h_PIDK_2.Fill(t.Lb_ISOLATION_PIDK)
                    if t.Lb_ISOLATION_PIDK2>4:
                        h_PIDK2_2.Fill(t.Lb_ISOLATION_PIDK2)


    c = r.TCanvas('c','c',1500,1500)
    c.Divide(4,2)
    c.cd(1)
    h_BDT.Draw()
    h_BDT_sel.SetFillStyle(3004)
    h_BDT_sel.SetFillColor(r.kRed)
    h_BDT_sel.Draw('same')
    print('Fraction of events with ISO_BDT > 0.35: ', h_BDT_sel.GetEntries()*1./h_BDT.GetEntries())
    c.cd(2)
    h_BDT2.Draw()
    h_BDT2_sel.SetFillStyle(3004)
    h_BDT2_sel.SetFillColor(r.kRed)
    h_BDT2_sel.Draw('same')
    print('Fraction of events with ISO_BDT2 > 0.2: ', h_BDT2_sel.GetEntries()*1./h_BDT2.GetEntries())
    c.cd(3)
    h_TrueD1.Draw()
    h_TrueD1_1.SetFillStyle(3004)
    h_TrueD1_1.SetFillColor(r.kRed)
    h_TrueD1_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2: ',  h_TrueD1_1.GetEntries()*1/h_TrueD1.GetEntries())
    c.cd(4)
    h_TrueD2.Draw()
    h_TrueD2_1.SetFillStyle(3004)
    h_TrueD2_1.SetFillColor(r.kRed)
    h_TrueD2_1.Draw('same')
    c.cd(5)
    h_PIDK.Draw()
    h_PIDK_1.SetFillStyle(3004)
    h_PIDK_1.SetFillColor(r.kRed)
    h_PIDK_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2 and PIDK>4: ', h_PIDK_2.GetEntries()*1./ h_PIDK_1.GetEntries())
    c.cd(6)
    h_PIDK2.Draw()
    h_PIDK2_1.SetFillStyle(3004)
    h_PIDK2_1.SetFillColor(r.kRed)
    h_PIDK2_1.Draw('same')
    print('Fraction of evts with ISO_BDT>0.35 and ISO_BDT2>0.2 and PIDK2>4: ', h_PIDK2_2.GetEntries()*1./ h_PIDK2_1.GetEntries())
    c.Update()
    c.Draw()

    return 


'''
GetTotalNumberMCevts()
GetTotalNumber2body()
GetTotalNumberMbody()
print()
GetTotalNumberMCevts_TMatch()
print()
GetTotalNumberMCevts_Trigger()
print()
GetTotalNumberMCevts_FullPreselection()

GetTotalNumberMCevts_TMatch_Iso()
GetTotalNumber2body_TMatch_Iso()
GetTotalNumberMbody_TMatch_Iso()
print()
GetTotalNumberMCevts_TMatch_Kenriched()
GetTotalNumber2body_TMatch_Kenriched()
GetTotalNumberMbody_TMatch_Kenriched()
print()
Check_mbody()
'''
