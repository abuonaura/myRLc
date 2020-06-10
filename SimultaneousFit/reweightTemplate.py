import formFactors
import ROOT as r

def theorDistr():
   hMuMean,hTauMean,hBoundries,c = formFactors.makeDistrWithParamVariation()
   return {"tau":{"max":hBoundries[0],"min":hBoundries[1],"mean":hTauMean}, "mu":{"max":hBoundries[2],"min":hBoundries[3],"mean":hMuMean}}

def mcDistr(tree,channel):
        histMC = r.TH1D("h_"+channel,"qem_"+channel,4,-2,14)
        q2= "(-((Lb_TRUEP_X-Lc_TRUEP_X)**2+(Lb_TRUEP_Y-Lc_TRUEP_Y)**2+(Lb_TRUEP_Z-Lc_TRUEP_Z)**2)+(Lb_TRUEP_E-Lc_TRUEP_E)**2)/1000000"
        tree.Draw(q2+">>h_"+channel,"Lc_BKGCAT<15&&Lb_BKGCAT<45")
        histMC = r.gPad.GetPrimitive("h_"+channel)
        histMC.Scale(1.0/histMC.Integral("width"))
        return histMC

def weightFitVar():
    outputFile = r.TFile.Open("DemoHistosLattice.root","RECREATE")
    templ = {"mu":"mu","tau":"tau"}
    histMC,hKi3 = {},{}
    hLattice = theorDistr()
    for iT in templ:
     f = r.TFile("/home/iaro/Work/LHCb/Analysis/dataMC1901/17Apr19/Lb_Lc"+iT+"nu_PID_reduced_preselected.root")
     t = f.Get("DecayTree")
     histMC[iT] = mcDistr(t,iT)
     hKi3[iT]={}
     for iL in hLattice[iT]:
        hKi3[iT][iL] = r.TH3D("h_w_"+iT+"_"+iL,"qem_w_"+iT+"_"+iL,4,-2,14,10,0,2600,10,-2,14)
        hLattice[iT][iL].Divide(histMC[iT])     
        for iE in t:
            if iE.Lc_BKGCAT>10 or iE.Lb_BKGCAT>40:continue
            coef = hLattice[iT][iL].GetBinContent(hLattice[iT][iL].FindBin(iE.FitVar_q2_mLc/1000000))
            hKi3[iT][iL].Fill(iE.FitVar_q2_mLc/1000000,iE.FitVar_El_mLc,iE.FitVar_Mmiss2_mLc/1000000,coef)
        #hKi3[iT][iL].Scale(1.0/hKi3[iT][iL].Integral("width"))
        hKi3[iT][iL].SetDirectory(0)
        outputFile.cd()
        hKi3[iT][iL].Write()
     f.Close()
    outputFile.Close()
    return hKi3
