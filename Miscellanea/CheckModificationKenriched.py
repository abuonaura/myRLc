import ROOT as r

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'

def SetNonNullBinContent(h):
    for i in range(h.GetNbinsX()):
        for j in range(h.GetNbinsY()):
            for k in range(h.GetNbinsZ()):
                if h.GetBinContent(i+1,j+1,k+1) <= 0:
                    h.SetBinContent(i+1,j+1,k+1,1E-6)
    return


if __name__== "__main__":
    dtype = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
    polarities=['MagUp','MagDown']
    
    histos = []
    histos_cut = []
    
    hTOT = r.TH3F('hTOT',"qem",4,-2,14,10,0,2600,10,-2,14)
    hTOT_cut = r.TH3F('hTOT_cut',"qem",4,-2,14,10,0,2600,10,-2,14)

    for dt in dtype:
        print(dt)
        h = r.TH3F('h_'+dt,"qem_"+dt,4,-2,14,10,0,2600,10,-2,14)
        h_cut = r.TH3F('h_cut_'+dt,"qem_"+dt,4,-2,14,10,0,2600,10,-2,14)
        h.SetDirectory(0)
        h_cut.SetDirectory(0)
        for polarity in polarities:
            inputFile = filedir+'Lb_'+dt+'_'+polarity+'_full_preselected_Kenr.root'
            f = r.TFile(inputFile,'READ')
            t = f.Get('DecayTree')
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                weight=t.Event_PIDCalibEffWeight
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    hTOT.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    if t.mu_ProbNNghost<0.2:
                        h_cut.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                        hTOT_cut.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)

        histos.append(h)    
        histos_cut.append(h_cut)    

    for h in histos:
        SetNonNullBinContent(h)
    for h in histos_cut:
        SetNonNullBinContent(h_cut)
    SetNonNullBinContent(hTOT)
    SetNonNullBinContent(hTOT_cut)

    #I believe that comparing the tot histos won't show a dramatical change because it does not take into account the different size of the samples
    hTOT_q2 = hTOT.ProjectionX()
    hTOT_El = hTOT.ProjectionY()
    hTOT_Mmiss2 = hTOT.ProjectionZ()
    hTOT_cut_q2 = hTOT_cut.ProjectionX()
    hTOT_cut_El = hTOT_cut.ProjectionY()
    hTOT_cut_Mmiss2 = hTOT_cut.ProjectionZ()
    c = r.TCanvas('c','c',1500,500)
    c.Divide(3,1)
    c.cd(1)
    hTOT_q2.SetLineColor(r.kBlue)
    hTOT_cut_q2.SetLineColor(r.kRed)
    hTOT_q2.Draw('hist')
    hTOT_cut_q2.Draw('hist same')
    c.cd(2)
    hTOT_El.SetLineColor(r.kBlue)
    hTOT_cut_El.SetLineColor(r.kRed)
    hTOT_El.Draw('hist')
    hTOT_cut_El.Draw('hist same')
    c.cd(3)
    hTOT_Mmiss2.SetLineColor(r.kBlue)
    hTOT_cut_Mmiss2.SetLineColor(r.kRed)
    hTOT_Mmiss2.Draw('hist')
    hTOT_cut_Mmiss2.Draw('hist same')

    #I now compare what happens to all the different samples (I do not expect major changes in lcmunu, lctaunu samples)
    canvas = []
    for i,h in enumerate(histos):
        h_cut = histos_cut[i]
        h_q2 = h.ProjectionX()
        h_El = h.ProjectionY()
        h_Mmiss2 = h.ProjectionZ()
        h_cut_q2 = h_cut.ProjectionX()
        h_cut_El = h_cut.ProjectionY()
        h_cut_Mmiss2 = h_cut.ProjectionZ()
        c1 = r.TCanvas('c1_'+str(i),'c1_'+str(i),1500,500)
        c1.Divide(3,1)
        c1.cd(1)
        h_q2.SetLineColor(r.kBlue)
        h_cut_q2.SetLineColor(r.kRed)
        h_q2.Draw('hist')
        h_cut_q2.Draw('hist same')
        c1.cd(2)
        h_El.SetLineColor(r.kBlue)
        h_cut_El.SetLineColor(r.kRed)
        h_El.Draw('hist')
        h_cut_El.Draw('hist same')
        c1.cd(3)
        h_Mmiss2.SetLineColor(r.kBlue)
        h_cut_Mmiss2.SetLineColor(r.kRed)
        h_Mmiss2.Draw('hist')
        h_cut_Mmiss2.Draw('hist same')
        c1.SaveAs('plots/FitVars_'+dtype[i]+'.png')
        canvas.append(c1)


