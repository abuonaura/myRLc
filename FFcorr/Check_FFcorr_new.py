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
    h.ProjectionX().Draw('hist sames')
    c.cd(2)
    h.ProjectionY().Draw('hist sames')
    c.cd(3)
    h.ProjectionZ().Draw('hist sames')
    return c

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

samples    = ['Lcmunu','Lctaunu','Lc2593munu','Lc2625munu','Lc2593taunu','Lc2625taunu']
polarities = ['MagUp','MagDown']
categories = ['Isolated','Kenriched','Lcpipi']


h_iso   = r.TH3F('h_iso',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_ffiso = r.TH3F('h_ffiso',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_Kenr   = r.TH3F('h_Kenr',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_ffKenr = r.TH3F('h_ffKenr',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_Lcpp   = r.TH3F('h_Lcpp',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
h_ffLcpp = r.TH3F('h_ffLcpp',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)

for sample in samples:
    print('Sample: ',sample)
    h_iso.Reset()
    h_ffiso.Reset()
    h_Kenr.Reset()
    h_ffKenr.Reset()
    h_Lcpp.Reset()
    h_ffLcpp.Reset()
    for polarity in polarities:
        print('  ',polarity)
        fname      = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full.root'
        tname      = 'tupleout/DecayTree'
        fname_psel = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_preselectionVars.root'
        tname_psel = 'DecayTree'
        fname_PID  = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_PIDCalib.root'
        tname_PID  = 'tupleout/DecayTree'
        fname_test = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/Lb_'+sample+'_'+polarity+'_full_FFcorrections.root'
        tname_test = 'DecayTree'


        vars2read = ['FinalSel','isIsolated','isKenriched','isLcpipi','Lb_True_Q2_mu','Lb_True_Costhetal_mu','TruthMatch','Event_FFcorr','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc','Event_PIDCalibEffWeight']

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

        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.FinalSel!=1: continue
            w = t.Event_PIDCalibEffWeight
            w_ff = t.Event_FFcorr
            if t.isIsolated==1:
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h_iso.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
                    h_ffiso.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w*w_ff)
            if t.isKenriched==1:
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h_Kenr.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
                    h_ffKenr.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w*w_ff)
            if t.isLcpipi==1:
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h_Lcpp.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
                    h_ffLcpp.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w*w_ff)

    h_iso = ScaleHisto(h_iso,1)
    h_ffiso = ScaleHisto(h_ffiso,1)

    h_Kenr = ScaleHisto(h_Kenr,1)
    h_ffKenr = ScaleHisto(h_ffKenr,1)

    h_Lcpp = ScaleHisto(h_Lcpp,1)
    h_ffLcpp = ScaleHisto(h_ffLcpp,1)

    c = PlotTemplates(h_iso,'c')
    c = PlotTemplatesSameCanvas(c,h_ffiso,r.kOrange-3)
    c.SaveAs('plots2/FitVars_'+sample+'_iso.png')
    
    c1 = PlotTemplates(h_Kenr,'c1')
    c1 = PlotTemplatesSameCanvas(c1,h_ffKenr,r.kOrange-3)
    c1.SaveAs('plots2/FitVars_'+sample+'_Kenr.png')
    
    c2 = PlotTemplates(h_Lcpp,'c2')
    c2 = PlotTemplatesSameCanvas(c2,h_ffLcpp,r.kOrange-3)
    c2.SaveAs('plots2/FitVars_'+sample+'_Lcpipi.png')
