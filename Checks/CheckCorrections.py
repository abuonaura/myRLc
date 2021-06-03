import ROOT as r

folder = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
mcsamples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds',
             'Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds']
polarities = ['MagUp','MagDown']

def PutTogetherPolarityHistos(h):
    h_new = {}
    for mcsample in mcsamples:
        h_new[mcsample] = h[mcsample]['MagUp']
        h_new[mcsample].Add(h[mcsample]['MagDown'])
        h_new[mcsample].SetDirectory(0)
    return h_new

def GetTemplateHistos(cut):
    h_q2 = {mcsample: {polarity: r.TH1F('h_q2_'+mcsample+'_'+polarity,'',40,-2,12) for polarity in polarities} for mcsample in mcsamples}
    h_MM2 = {mcsample: {polarity: r.TH1F('h_Mmiss2_'+mcsample+'_'+polarity,'',100,-2,12) for polarity in polarities} for mcsample in mcsamples}
    h_El = {mcsample: {polarity: r.TH1F('h_El2_'+mcsample+'_'+polarity,'',100,0,2600) for polarity in polarities} for mcsample in mcsamples}
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            fPID = r.TFile(folder+mcsample+'_'+polarity+'_full_PIDCalib.root','READ')
            tPID = fPID.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            startt.AddFriend(tPID)
            name = mcsample+'_'+polarity
            startt.Draw('FitVar_q2_mLc/1.E6>>q2_'+name+'(40,-2,14)','w_LbCorr*Event_PIDCalibEffWeight*(FitVar_q2_mLc/1.E6>-2 && FitVar_q2_mLc/1.E6<14'+cut+')')
            h_q2[mcsample][polarity] = r.gPad.GetPrimitive('q2_'+name) 
            startt.Draw('FitVar_Mmiss2_mLc/1.E6>>MM2_'+name+'(100,-2,14)','w_LbCorr*Event_PIDCalibEffWeight*(FitVar_Mmiss2_mLc/1.E6>-2 && FitVar_Mmiss2_mLc/1.E6<14 '+cut+')')
            h_MM2[mcsample][polarity] = r.gPad.GetPrimitive('MM2_'+name) 
            startt.Draw('FitVar_El_mLc/1.E6>>El_'+name+'(100,0,2600)','w_LbCorr*Event_PIDCalibEffWeight*(FitVar_El_mLc>0 && FitVar_El_mLc<2600 ' +cut+')')
            h_El[mcsample][polarity] = r.gPad.GetPrimitive('El_'+name) 
            h_q2[mcsample][polarity].SetDirectory(0)
            h_MM2[mcsample][polarity].SetDirectory(0)
            h_El[mcsample][polarity].SetDirectory(0)

    h_q2_1 = PutTogetherPolarityHistos(h_q2)
    h_MM2_1 = PutTogetherPolarityHistos(h_MM2)
    h_El_1 = PutTogetherPolarityHistos(h_El)
    return h_q2_1,h_MM2_1, h_El_1

def PlotSuperimposed(h,h_TM,h_Tr,h_Presel,h_Final):
    c = r.TCanvas()
    h.SetLineColor(r.Black)
    h_TM.SetLineColor(r.kOrange+2)
    h_Tr.SetLineColor(r.kGreen-3)
    h_Presel.SetLineColor(r.kAzure+9)
    h_Final.SetLineColor(r.kRed)
    h.Draw('hist')
    h_TM.Draw('hist same')
    h_Tr.Draw('hist same')
    h_Presel.Draw('hist same')
    h_Final.Draw('hist same')
    l = r.TLegend()
    l.AddEntry(h,'start','l')
    l.AddEntry(h_TM,'TruthMatched','l')
    l.AddEntry(h_Trigger,'Trigger','l')
    l.AddEntry(h_Presel,'Presel','l')
    l.AddEntry(h_Final,'Final Presel','l')
    l.Draw('same')
    return c

h_q2,h_MM2, h_El = GetTemplateHistos('')
h_q2_TM,h_MM2_TM, h_El_TM = GetTemplateHistos('&&TruthMatch==1')
h_q2_Tr,h_MM2_Tr, h_El_Tr = GetTemplateHistos('&&TruthMatch==1&&Trigger==1&&MoreHLT2Sel==1')
h_q2_Presel,h_MM2_Presel, h_El_Presel = GetTemplateHistos('&&TruthMatch==1&&Trigger==1&&MoreHLT2Sel==1&&PIDCalib==1&&LcMass==1')
h_q2_Final,h_MM2_Final, h_El_Final = GetTemplateHistos('&&FinalSel==1')
h_q2_iso,h_MM2_iso, h_El_iso = GetTemplateHistos('&&FinalSel==1&&isIsolated==1')
h_q2_Kenr,h_MM2_Kenr, h_El_Kenr = GetTemplateHistos('&&FinalSel==1&&isKenriched==1')
h_q2_Lcpipi,h_MM2_Lcpipi, h_El_Lcpipi = GetTemplateHistos('&&FinalSel==1&&isLcpipi==1')


