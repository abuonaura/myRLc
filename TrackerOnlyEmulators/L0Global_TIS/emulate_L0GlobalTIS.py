import ROOT as r
import math as m
from array import array
import sys

def main(MC_fileName = "/home/iaro/Work/LHCb/Analysis/share_scripts/fromSimone/L0HLT1emulationVars/Lb_Lcmunu_MagUp.root"):
    
    JpsiKFile = {"file":"L0TIS_efficiency_2D_2016Data.root",
                "tree":"TISTOS"
                }
    reweighter_file = r.TFile(JpsiKFile["file"])
    reweighter_histo = reweighter_file.Get(JpsiKFile["tree"])

    fin = r.TFile.Open(MC_fileName)
    tin = fin.Get('tupleout/DecayTree')
    tin.SetBranchStatus('*',0)
    tin.SetBranchStatus('Lb_TRUEP*',1)
    fout = r.TFile.Open(MC_fileName[:-5]+'_wL0TISEmulation.root','RECREATE')
    tout = r.TTree("DecayTree", "DecayTree")

    Lb_L0Global_TIS_emulated = array('i',[0])
    br_L0TIS_weight = tout.Branch('Lb_L0Global_TIS_emulated', Lb_L0Global_TIS_emulated, 'Lb_L0Global_TIS_emulated[1]/I')

    rndm = r.TRandom3(1)
    print("Processing events ...")
    for ev in tin:
            s = rndm.Uniform(0,1.)
            if ev.Lb_TRUEP_Z>0:
                    xbin = reweighter_histo.GetXaxis().FindBin(m.log(ev.Lb_TRUEP_Z))
                    ybin = reweighter_histo.GetYaxis().FindBin(m.log(ev.Lb_TRUEPT))
                    weight = reweighter_histo.GetBinContent(xbin,ybin)
                    #print (weight,s)
                    if weight>s:
                        Lb_L0Global_TIS_emulated[0] = 1
                    else:
                        Lb_L0Global_TIS_emulated[0] = 0
            else: Lb_L0Global_TIS_emulated[0] = 0
            tout.Fill()
    print("All events processed.")

    tout.Write()
    fout.Close()
    fin.Close()
