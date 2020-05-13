import ROOT as r
import math as m
from array import array

def PtDes(T):
    PT=[]
    for iND in {'p','pi','K'}:
     PT.append(T.GetLeaf(iND+"_PT").GetValue())
    PT.sort()
    if PT[2]>600 and PT[1]>600:
        return r.kTRUE

def PDes(T):
    P=[]
    for iND in {'p','pi','K'}:
     P.append(T.GetLeaf(iND+"_P").GetValue())
    P.sort()
    if P[2]>5000 and P[1]>5000:
        return r.kTRUE
    
def SUMPTDes(T):
 if T.nVeloClusters >= 6000 or T.nITClusters >=3000 or T.nOTClusters>=15000 or T.nVeloClusters <= 50 or T.nITClusters <=50 or T.nOTClusters<=50:
     return r.kFALSE
 decision = 0
 SUMPT=[]
 name=[]
 for iND in {'p','K'}:
     for iD in {'K','pi'}:
         if iND!=iD:
            PX = T.GetLeaf(iND+"_PX").GetValue() +T.GetLeaf(iD+"_PX").GetValue()
            PY = T.GetLeaf(iND+"_PY").GetValue() +T.GetLeaf(iD+"_PY").GetValue()
            SUMPT.append(m.sqrt(PX*PX+PY*PY))
            name.append([iND,iD,m.sqrt(PX*PX+PY*PY)])
 SUMPT.sort()
 naming={'p':'1','K':'2','pi':'3'}
 if SUMPT[2]>2000:
     for iNr in name:
         if iNr[2]==SUMPT[2]:
            if T.GetLeaf(iNr[0]+"_P").GetValue()>5000 and  T.GetLeaf(iNr[1]+"_P").GetValue()>5000:
                if T.GetLeaf(iNr[0]+"_PT").GetValue()>600 and  T.GetLeaf(iNr[1]+"_PT").GetValue()>600:
                    ######################################
                    if T.GetLeaf("Lb_IPCHI2_OWNPV_COMB_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>4:
                         
                        if T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[0]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[0]]).GetValue()< 2.5 and T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[1]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[1]]).GetValue()< 2.5:
                            if T.GetLeaf("Lb_HLt1TwoTrackMVAEmulations_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>0.95:
                                if T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[0]]).GetValue()<0.2 and T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[1]]).GetValue()<0.2:
                                    if T.GetLeaf(iNr[0]+"_TRACK_nTTHits").GetValue()>2 and T.GetLeaf(iNr[1]+"_TRACK_nTTHits").GetValue()>2:
                                        decision = 1
                    ########################################3
 if decision==0:
     if SUMPT[1]>2000:
      for iNr in name:
         if iNr[2]==SUMPT[1]:
            if T.GetLeaf(iNr[0]+"_P").GetValue()>5000 and  T.GetLeaf(iNr[1]+"_P").GetValue()>5000:
                if T.GetLeaf(iNr[0]+"_PT").GetValue()>600 and  T.GetLeaf(iNr[1]+"_PT").GetValue()>600:
                    ######################################
                    if T.GetLeaf("Lb_IPCHI2_OWNPV_COMB_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>4:
                         
                        if T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[0]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[0]]).GetValue()< 2.5 and T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[1]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[1]]).GetValue()< 2.5:
                            if T.GetLeaf("Lb_HLt1TwoTrackMVAEmulations_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>0.95:
                                if T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[0]]).GetValue()<0.2 and T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[1]]).GetValue()<0.2:
                                    if T.GetLeaf(iNr[0]+"_TRACK_nTTHits").GetValue()>2 and T.GetLeaf(iNr[1]+"_TRACK_nTTHits").GetValue()>2:
                                        decision = 1
                                
 if decision==0:
     if SUMPT[0]>2000:
      for iNr in name:
         if iNr[2]==SUMPT[0]:
            if T.GetLeaf(iNr[0]+"_P").GetValue()>5000 and  T.GetLeaf(iNr[1]+"_P").GetValue()>5000:
                if T.GetLeaf(iNr[0]+"_PT").GetValue()>600 and  T.GetLeaf(iNr[1]+"_PT").GetValue()>600:
                    ######################################
                    if T.GetLeaf("Lb_IPCHI2_OWNPV_COMB_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>4:
                         
                        if T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[0]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[0]]).GetValue()< 2.5 and T.GetLeaf("Lb_TRACK_CHI2_DAU_"+naming[iNr[1]]).GetValue()/T.GetLeaf("Lb_TRACK_NDOF_DAU_"+naming[iNr[1]]).GetValue()< 2.5:
                           if T.GetLeaf("Lb_HLt1TwoTrackMVAEmulations_"+naming[iNr[0]]+"_"+naming[iNr[1]]).GetValue()>0.95:
                                if T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[0]]).GetValue()<0.2 and T.GetLeaf("Lb_TRACK_GHOSTPROB_DAU_"+naming[iNr[1]]).GetValue()<0.2:
                                    if T.GetLeaf(iNr[0]+"_TRACK_nTTHits").GetValue()>2 and T.GetLeaf(iNr[1]+"_TRACK_nTTHits").GetValue()>2:
                                        decision = 1                               
                                
                                
 if decision==1: return r.kTRUE

def main(MC_fileName = "Lb_Lcmunu_MagUp.root"):
    F=r.TFile(MC_fileName)
    T=F.Get("tupleout/DecayTree") 
    fHLT1TwoTracks = r.TFile(MC_fileName[:-5]+"_wHLT1TwoTracksEmulation.root","RECREATE")   
    tHLT1TwoTracks = r.TTree("DecayTree","HLT1 Two Tracks Emulations")
    HLT1TwoTracksDecision = array('i',[0])
    tHLT1TwoTracks.Branch('Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS',HLT1TwoTracksDecision,'HLT1TwoTracksDecision[1]/I')

    hP=r.TH1F("Lc_P","Lc_p",100,0,1000000)
    for iE in range(T.GetEntries()):
        T.GetEntry(iE)
        HLT1TwoTracksDecision[0] = 0
        if T.GetLeaf("nSPDHits").GetValue()>450:
            tHLT1TwoTracks.Fill()
            continue
        if T.GetLeaf("Lc_L0HadronDecision_TOS").GetValue()>0 or T.GetLeaf("Lb_L0Global_TIS").GetValue()>0:
            if PtDes(T)==r.kTRUE and PDes(T)==r.kTRUE and SUMPTDes(T)==r.kTRUE:
                hP.Fill(T.Lc_P)
                HLT1TwoTracksDecision[0] = 1
        tHLT1TwoTracks.Fill()
    tHLT1TwoTracks.Write()
    fHLT1TwoTracks.Write()
    fHLT1TwoTracks.Close()
 
