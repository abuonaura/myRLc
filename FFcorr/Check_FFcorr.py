import ROOT as r
import uproot4

def PlotHisto(h,name):
    c = r.TCanvas(name,'',500,500)
    h.SetLineColor(r.kAzure-3)
    h.Draw('hist')
    return c

def PlotHisto2D(h,name):
    c = r.TCanvas(name,'',500,500)
    h.Draw('colz')
    return c

def PlotHistoSameCanvas(c,h,color):
    h.SetLineColor(color)
    h.Draw('hist same')
    return c

def PlotTemplates(h,name):
    h.SetLineColor(r.kAzure-3)
    c = r.TCanvas(name,'',1500,500)
    c.Divide(3,1)
    c.cd(1)
    h.ProjectionX().Draw('hist')
    c.cd(2)
    h.ProjectionY().Draw('hist')
    c.cd(3)
    h.ProjectionZ().Draw('hist')
    return c

def PlotTemplatesSameCanvas(c,h,color):
    h.SetLineColor(color)
    c.cd(1)
    h.ProjectionX().Draw('hist same')
    c.cd(2)
    h.ProjectionY().Draw('hist same')
    c.cd(3)
    h.ProjectionZ().Draw('hist same')
    return c

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

samples     = ['Lcmunu']
#samples    = ['Lcmunu','Lctaunu','Lc2593munu','Lc2625munu','Lc2593taunu','Lc2625taunu']
#samples    = ['Lctaunu','Lc2593munu','Lc2625munu','Lc2593taunu','Lc2625taunu']
polarities = ['MagUp','MagDown']
categories = ['Isolated','Kenriched','Lcpipi']

h_q2_tmatched = r.TH1F('h_q2_tmatched','q^{2} for TMatched sample',30,-2E6,14E6)
h_q2_fsel     = r.TH1F('h_q2_fsel','q^{2} for sample passing full preselection',30,-2E6,14E6)
h_q2_iso      = r.TH1F('h_q2_iso','q^{2} for sample passing isolation cut',30,-2E6,14E6)

h_q2_cthl_iso = r.TH2F('h_q2_cthl_iso','q^{2} and cos #theta_{l} for isolated sample',30,-2E6,14E6,50,-1.5,1.5)
h_q2_ff_iso = r.TH2F('h_q2_ff_iso','q^{2} and ff correction for iso sample',30,-2E6,14E6,50,0,1.5)
h_ff_cthl_iso = r.TH2F('h_ff_cthl_iso','cos #theta_{l} and ff correction for iso sample',50,-1.5,1.5,50,0,1.5)

h_cthl_tmatched = r.TH1F('h_cthl_tmatched','cos #theta_{l} for TMatched sample',50,-2,2)
h_cthl_fsel     = r.TH1F('h_cthl_fsel','cos #theta_{l} for sample passing full preselection',50,-2,2)
h_cthl_iso      = r.TH1F('h_cthl_iso','cos #theta_{l} for sample passing isolation cut',50,-2,2)

h_templ_iso    = r.TH3F('h_templ_iso',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_templ_ffcorr = r.TH3F('h_templ_ffcorr',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)

for sample in samples:
    print(sample)
    h_q2_tmatched.Reset()
    h_q2_fsel.Reset()
    h_q2_iso.Reset()
    h_q2_cthl_iso.Reset()
    h_q2_ff_iso.Reset()
    h_ff_cthl_iso.Reset()
    h_cthl_tmatched.Reset()
    h_cthl_fsel.Reset()
    h_cthl_iso.Reset()
    h_templ_iso.Reset()
    h_templ_ffcorr.Reset()

    for polarity in polarities:
        print(polarity)
        print(h_q2_fsel)
        fname      = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full.root'
        tname      = 'tupleout/DecayTree'
        fname_psel = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_preselectionVars.root'
        tname_psel = 'DecayTree'
        fname_PID  = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_PIDCalib.root'
        tname_PID  = 'tupleout/DecayTree'
        fname_test = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_FFcorrections.root'
        tname_test = 'DecayTree'


        vars2read = ['FinalSel','isIsolated','Lb_True_Q2_mu','Lb_True_Costhetal_mu','TruthMatch','Event_FFcorr','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc','Event_PIDCalibEffWeight']

        f      = r.TFile(fname, 'READ')
        f_psel = r.TFile(fname_psel,'READ')
        f_PID  = r.TFile(fname_PID,'READ')
        f_test = r.TFile(fname_test,'READ')

        t      = f.Get(tname)
        t_psel = f_psel.Get(tname_psel)
        t_PID  = f_PID.Get(tname_PID)
        t_test = f_test.Get(tname_test)

        t.AddFriend(t_psel)
        t.AddFriend(t_PID)
        t.AddFriend(t_test)

        #Switch off all branches to speed up
        t.SetBranchStatus('*',0)
        #and turn on only those needed
        for b in vars2read:
            t.SetBranchStatus(b,1)



        #for i in range(t.GetEntries()):
        for i in range(int(t.GetEntries()/2)):
            t.GetEntry(i)
            if t.FinalSel != 1:
                continue
            else:
                w = t.Event_PIDCalibEffWeight
                w_ff = t.Event_FFcorr
                h_q2_fsel.Fill(t.Lb_True_Q2_mu,w)
                h_cthl_fsel.Fill(t.Lb_True_Costhetal_mu,w)
                if t.isIsolated==1:
                    h_q2_iso.Fill(t.Lb_True_Q2_mu,w)
                    h_cthl_iso.Fill(t.Lb_True_Costhetal_mu,w)
                    h_q2_cthl_iso.Fill(t.Lb_True_Q2_mu,t.Lb_True_Costhetal_mu,w)
                    h_q2_ff_iso.Fill(t.Lb_True_Q2_mu,w_ff,w)
                    h_ff_cthl_iso.Fill(t.Lb_True_Costhetal_mu,w_ff,w)
                    h_templ_iso.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
                    h_templ_ffcorr.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w*w_ff)



    #c = PlotHisto(h_q2_tmatched,'c')
    c = PlotHisto(h_q2_fsel,'c')
    c = PlotHistoSameCanvas(c, h_q2_iso,r.kGreen+2)
    c.SaveAs('plots/q2_'+sample+'.png')

    #c1 = PlotHisto(h_cthl_tmatched,'c1')
    #c1 = PlotHistoSameCanvas(c1, h_cthl_fsel,r.kOrange-3)
    c1 = PlotHisto(h_cthl_fsel,'c1')
    c1 = PlotHistoSameCanvas(c1, h_cthl_iso,r.kGreen+2)
    c1.SaveAs('plots/cthl_'+sample+'.png')

    c2 = PlotTemplates(h_templ_iso,'c2')
    c2 = PlotTemplatesSameCanvas(c2,h_templ_ffcorr,r.kOrange-3)
    c2.SaveAs('plots/FitVars_'+sample+'.png')

    c3 = PlotHisto2D(h_q2_cthl_iso,'c3') 
    c3.SaveAs('plots/q2VScthl_'+sample+'.png')
    c4 = PlotHisto2D(h_q2_ff_iso,'c4') 
    c4.SaveAs('plots/q2VSFF_'+sample+'.png')
    c5 = PlotHisto2D(h_ff_cthl_iso,'c5') 
    c5.SaveAs('plots/cthlVSFF_'+sample+'.png')

    h_templ_iso_norm    = ScaleHisto(h_templ_iso,1)
    h_templ_ffcorr_norm = ScaleHisto(h_templ_ffcorr,1)
    c6 = PlotTemplates(h_templ_iso_norm,'c6')
    c6 = PlotTemplatesSameCanvas(c6,h_templ_ffcorr_norm,r.kOrange-3)
    c6.SaveAs('plots/FitVars_'+sample+'_norm.png')
