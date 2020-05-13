import ROOT as r
import os
from AddSweights import *

filedir = '/disk/lhcb_data2/RLcMuonic2016/'

def CalculateSelectionEfficiencies(dtype, polarity):
    fsel = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Selections.root','READ')
    tsel = fsel.Get('DecayTree')

    nTriggerOK, nLcMOK, nMuOK, nPIDOK, nTOTpresOK, nBDTOK, nFinalOK, nISO, nKenr =0,0,0,0,0,0,0,0,0
    effTrigger, effLcM, effPID, effTOTpres, effBDT, effFinal, effISO, effKenr =0,0,0,0,0,0,0,0
    #for i in range(0,500000):
    for i in range(tsel.GetEntries()):
        tsel.GetEntry(i)
        if tsel.Trigger_Final==1:
            nTriggerOK+=1
            if tsel.LcMass==1:
                nLcMOK+=1
                if tsel.PIDcut==1:
                    nMuOK+=1
                    if tsel.PIDCalibCut==1:
                        nPIDOK+=1
                        if tsel.TOTpreselection==1:
                            nTOTpresOK+=1
                            if tsel.PassBDT==1:
                                nBDTOK+=1
                                if tsel.FinalSelection==1:
                                    nFinalOK+=1
                                    if tsel.isIsolated==1 and tsel.isKenriched==0:
                                        nISO+=1
                                    if tsel.isKenriched==1 and tsel.isIsolated==0:
                                        nKenr+=1

    print(nTriggerOK,nLcMOK,nMuOK,nPIDOK,nTOTpresOK,nBDTOK)
    effTrigger = nTriggerOK/tsel.GetEntries()
    effLcM = nLcMOK/nTriggerOK
    effPID = nPIDOK/nLcMOK
    effTOTpres = nTOTpresOK/nPIDOK
    effBDT = nBDTOK/nTOTpresOK
    effFinal = nFinalOK/nBDTOK
    effISO = nISO/nFinalOK
    effKenr = nKenr/nFinalOK
    print(nISO, nFinalOK, nKenr)
    print('Trigger efficiency: {:.3f}'.format(effTrigger))
    print('LcMass efficiency: {:.3f}'.format(effLcM))
    print('Cuts on PID vars efficiency: {:.3f}'.format(effPID))
    print('Total preselection efficiency: {:.3f}'.format(effTOTpres))
    print('BDT efficiency: {:.3f}'.format(effBDT))
    print('Final selection efficiency: {:.3f}'.format(effFinal))
    print('Isolation efficiency: {:.3f}'.format(effISO))
    print('Kenriched efficiency: {:.3f}'.format(effKenr))

def CalculateSelectionEfficienciesDf(dtype, polarity):
    fsel = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Df_preselection.root','READ')
    tsel = fsel.Get('DecayTree')

    nTriggerOK, nLcMOK, nMuOK, nPIDOK, nTOTpresOK, nBDTOK, nFinalOK, nISO, nKenr =0,0,0,0,0,0,0,0,0
    effTrigger, effLcM, effPID, effTOTpres, effBDT, effFinal, effISO, effKenr =0,0,0,0,0,0,0,0
    for i in range(tsel.GetEntries()):
    #for i in range(0,500000):
        tsel.GetEntry(i)
        if tsel.Trigger==1:
            nTriggerOK+=1
            if tsel.LcMass==1:
                nLcMOK+=1
                if tsel.MuCuts==1:
                    nMuOK+=1
                    if tsel.PIDCalib==1:
                        nPIDOK+=1
                        if tsel.Preselection==1:
                            nTOTpresOK+=1
                            if tsel.PassBDT==1:
                                nBDTOK+=1
                                if tsel.FinalSel==1:
                                    nFinalOK+=1
                                    if tsel.isIsolated==1 and tsel.isKenriched==0:
                                        nISO+=1
                                    if tsel.isKenriched==1 and tsel.isIsolated==0:
                                        nKenr+=1
    print(nTriggerOK,nLcMOK,nMuOK,nPIDOK,nTOTpresOK,nBDTOK)
    effTrigger = nTriggerOK/tsel.GetEntries()
    effLcM = nLcMOK/nTriggerOK
    effPID = nPIDOK/nLcMOK
    effTOTpres = nTOTpresOK/nPIDOK
    effBDT = nBDTOK/nTOTpresOK
    effFinal = nFinalOK/nBDTOK
    effISO = nISO/nFinalOK
    effKenr = nKenr/nFinalOK
    print(nISO, nFinalOK, nKenr)
    print('Trigger efficiency: {:.3f}'.format(effTrigger))
    print('LcMass efficiency: {:.3f}'.format(effLcM))
    print('Cuts on PID vars efficiency: {:.3f}'.format(effPID))
    print('Total preselection efficiency: {:.3f}'.format(effTOTpres))
    print('BDT efficiency: {:.3f}'.format(effBDT))
    print('Final selection efficiency: {:.3f}'.format(effFinal))
    print('Isolation efficiency: {:.3f}'.format(effISO))
    print('Kenriched efficiency: {:.3f}'.format(effKenr))


suffix = {'full':'.root', 'iso':'_iso.root','Kenriched':'_Kenr.root'}

def CreateFinalTree(dtype,polarity,sample):
    f = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'.root','READ')
    t = f.Get('tupleout/DecayTree')
    fsel = r.TFile(filedir+'Data/Lb_'+dtype+'_'+polarity+'_Df_preselection.root','READ')
    tsel = fsel.Get('DecayTree')
    t.AddFriend(tsel)
    ofname = filedir+'Data/Lb_'+dtype+'_'+polarity+'_preselected'+suffix[sample]
    of = r.TFile(ofname,'RECREATE')
    print('... Copying the tree ')
    if sample=='full':
        ot = t.CopyTree("FinalSel==1")
    if sample=='iso':
        ot = t.CopyTree("FinalSel==1&&isIsolated==1&&isKenriched==0")
    if sample=='Kenriched':
        ot = t.CopyTree("FinalSel==1&&isIsolated==0&&isKenriched==1")
    of.Write()
    print('Tree copied, moving to sweights')
    of.Close()
    swfname = ComputeSweights(ofname,dtype,polarity,sample)
    tname = 'DecayTree'
    variable = 'Lc_M'
    print('...Plotting')
    splotVariable(swfname, tname, dtype, polarity,sample,variable, 50, 2230, 2330, 'm_{#Lambda_{c}}')

    return


if __name__ == "__main__":
    #dtype = ['Data','DataSS','FakeMu','FakeMuSS']
    dtype = ['FakeMu','FakeMuSS']
    #polarities=['MagUp','MagDown']
    polarities=['MagDown']
    samples = ['full','iso','Kenriched']
    #samples = ['Kenriched']
    for polarity in polarities:
        print(polarity)
        for dt in dtype:
            print('   ', dt)
            #print(' >>>> Standard ROOT: ')
            #CalculateSelectionEfficiencies(dt,polarity)
            #print(' >>>> DataFrames: ')
            #CalculateSelectionEfficienciesDf(dt,polarity)
            for sample in samples:
                print('      ',sample)
                CreateFinalTree(dt,polarity,sample)

