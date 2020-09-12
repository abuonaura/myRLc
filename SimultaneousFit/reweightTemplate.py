import formFactors
import ROOT as r
import argparse

def theorDistr():
   hMuMean,hTauMean,hBoundries,c = formFactors.makeDistrWithParamVariation()
   return {"tau":{"max":hBoundries[0],"min":hBoundries[1],"mean":hTauMean}, "mu":{"max":hBoundries[2],"min":hBoundries[3],"mean":hMuMean}}

def mcDistr(channel):
    filedirMC = '/disk/lhcb_data2/RLcMuonic2016/GeneratorLevelMC/'
    if channel=="mu":
        fname = 'LcMuNu_gen.root'
    if channel=="tau":
        fname = 'LcTauNu_gen.root'
    f = r.TFile(filedirMC+fname,'READ')
    t = f.Get('MCDecayTreeTuple/MCDecayTree')
    print(t)
    histMC = r.TH1D("h_"+channel,"qem_"+channel,4,-2,14)
    q2= "(-((Lambda_b0_TRUEP_X-Lambda_cplus_TRUEP_X)**2+(Lambda_b0_TRUEP_Y-Lambda_cplus_TRUEP_Y)**2+(Lambda_b0_TRUEP_Z-Lambda_cplus_TRUEP_Z)**2)+(Lambda_b0_TRUEP_E-Lambda_cplus_TRUEP_E)**2)/1000000"
    #tree.Draw(q2+">>h_"+channel,"Lc_BKGCAT<15&&Lb_BKGCAT<45")
    c = r.TCanvas()
    t.Draw(q2+">>h_"+channel)
    c.SaveAs("histMC_"+channel+".png")
    histMC = r.gPad.GetPrimitive("h_"+channel)
    histMC.SetDirectory(0)
    return histMC

def weightFitVar(category,MCtype):
    outputFile = r.TFile.Open("RootFiles/DemoHistosLattice_"+category+"_"+MCtype+".root","RECREATE")
    if MCtype=='MCfull':
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
    else:
        filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/'
    if category=='Isolated':
        cat = 'iso'
    else:
        cat = 'Kenr'
    templ = {"mu":"mu","tau":"tau"}
    histMC,hKi3,weights = {},{},{}
    polarities=['MagUp','MagDown']
    hLattice = theorDistr()
    for hT in hLattice:
        print('hT: ', hT)
        for h in hLattice[hT]:
            print('h: ', h)
            hLattice[hT][h].Scale(1.0/ hLattice[hT][h].Integral("width"))
            print(hLattice[hT][h])
    for iT in templ:
        print('iT: ', iT)
        hKi3[iT]={}
        weights[iT]={}
        histMC[iT] = mcDistr(iT)
        histMC[iT].Scale(1.0/histMC[iT].Integral("width"))
        for iL in hLattice[iT]:
            hKi3[iT][iL] = r.TH3D("h_w_"+iT+"_"+iL,"qem_w_"+iT+"_"+iL,4,-2,14,10,0,2600,10,-2,14)
            weights[iT][iL]=[]
            print (iL)
            for polarity in polarities:
                if MCtype=='MCfull':
                    f = r.TFile('/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_Lc'+iT+'nu_'+polarity+'_full_preselected_'+cat+'_LbCorr.root')
                else:
                    f = r.TFile('/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/Lb_Lc'+iT+'nu_'+polarity+'_preselected_'+cat+'_LbCorr.root')
                t = f.Get("DecayTree")
                print(hLattice[iT][iL].GetEntries())
                hLattice1 = r.TH1D('hLattice1','',4,-2,14)
                hLattice1.Divide(hLattice[iT][iL],histMC[iT])     
                hLattice1.Draw()
                nentries=0
                for iE in t:
                    nentries+=1
                    coef = hLattice1.GetBinContent(hLattice1.FindBin(iE.FitVar_q2_mLc/1000000))
                    #print(coef)
                    weights[iT][iL].append(coef)
                    #coef=1
                    weight = iE.Event_PIDCalibEffWeight*iE.w_LbCorr
                    hKi3[iT][iL].Fill(iE.FitVar_q2_mLc/1000000,iE.FitVar_El_mLc,iE.FitVar_Mmiss2_mLc/1000000,coef*weight)
                f.Close()
                print('Nentries for polarity ' + polarity+ ': '+str(nentries))
            #hKi3[iT][iL].Scale(1.0/hKi3[iT][iL].Integral("width"))
            hKi3[iT][iL].SetDirectory(0)
            outputFile.cd()
            hKi3[iT][iL].Write()
    outputFile.Close()
    return hKi3, histMC, weights

def init():
    ap = argparse.ArgumentParser(description='Add Lb production corrections')
    ap.add_argument('-c','--category', type=str, dest='category', default=None)
    ap.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    ap.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    args = ap.parse_args()
    return args

def DrawLatticeHistos():
    hLattice = theorDistr()
    c = r.TCanvas('c','c')
    hLattice['tau']['max'].SetLineColor(r.kRed)
    hLattice['tau']['max'].SetFillColor(r.kWhite)
    hLattice['tau']['min'].SetLineColor(r.kBlue)
    hLattice['tau']['mean'].SetLineColor(r.kBlack)
    hLattice['tau']['max'].Draw()
    hLattice['tau']['mean'].Draw("sames")
    hLattice['tau']['min'].Draw("sames")
    leg = r.TLegend(0.2,0.4,0.4,0.6)
    leg.AddEntry(hLattice['tau']['max'],"max","L")
    leg.AddEntry(hLattice['tau']['mean'],"mean","L")
    leg.AddEntry(hLattice['tau']['min'],"min","L")
    leg.Draw()
    c.SaveAs('hLatticeTau.png')
    c1 = r.TCanvas('c1','c1')
    hLattice['mu']['max'].SetLineColor(r.kRed)
    hLattice['mu']['max'].SetFillColor(r.kWhite)
    hLattice['mu']['min'].SetLineColor(r.kBlue)
    hLattice['mu']['mean'].SetLineColor(r.kBlack)
    hLattice['mu']['max'].Draw()
    hLattice['mu']['mean'].Draw("sames")
    hLattice['mu']['min'].Draw("sames")
    leg1 = r.TLegend(0.2,0.4,0.4,0.6)
    leg1.AddEntry(hLattice['mu']['max'],"max","L")
    leg1.AddEntry(hLattice['mu']['mean'],"mean","L")
    leg1.AddEntry(hLattice['mu']['min'],"min","L")
    leg1.Draw()
    c1.SaveAs('hLatticeMu.png')
    c1.SaveAs('hLatticeMu.root')

if __name__== "__main__":
    args = init()
    category = args.category
    if args.MCfull==True:
        MCtype = 'MCfull'
    if args.MCTO==True:
        MCtype = 'MCTrackerOnly'

    print('Category: ', category)
    print('MCtype: ', MCtype)
    hKi3, histMC, weights=weightFitVar(category, MCtype)
    #DrawLatticeHistos()
