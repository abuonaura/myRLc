import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas, TMath, TLegend
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt

folderD = '/disk/lhcb_data2/RLcMuonic2016/Data/'
folderMISID = '/disk/lhcb_data2/RLcMuonic2016/MISID/OppositeSign/'
folderComb = '/disk/lhcb_data2/RLcMuonic2016/CombinatorialBkg/'
folderMC = '/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/'

datasamples = ['Data','MISID','Combinatorial']
mcsamples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds']
polarities = ['MagUp','MagDown']

#Variables to change to include/exclude some corrections
suffixl = []
weighComb = False
FFGstate  = False
FFEstateL = False
FFEstateH = True

def SetNonNullBinContent(h):
    for i in range(h.GetNbinsX()):
        for j in range(h.GetNbinsY()):
            for k in range(h.GetNbinsZ()):
                if h.GetBinContent(i+1,j+1,k+1) <= 0:
                    h.SetBinContent(i+1,j+1,k+1,1E-6)
    return


def Weight_MultibodyCharm_linear(w_mbody,alpha1):
    weight = 1+2*alpha1*w_mbody
    return weight

def Weight_MultibodyCharm_quadratic(w_mbody,alpha2):
    weight = (1-alpha2)+8*alpha2*w_mbody*w_mbody
    return weight

def Weight_Combinatorial(tree,HighMassCorrection):
    if HighMassCorrection == False:
        weight = getattr(tree,'w_MISID')*getattr(tree,'sw_sig')
    else:
        weight = getattr(tree,'w_MISID')*getattr(tree,'sw_sig')*getattr(tree,'w_comb')
    return weight

def Weight_MISID(tree):
    weight = getattr(tree,'w_recomu_CF')*10*getattr(tree,'sw_sig')
    return weight

def Weight_Data(tree):
    return getattr(tree,'sw_sig')

def GetWeight(sample,tree,MCtype):
    FF_Lbmunu      = FFGstate 
    FF_Lbtaunu     = FFGstate
    FF_Lb2625munu  = FFEstateH 
    FF_Lb2625taunu = FFEstateH 
    FF_Lb2593munu  = FFEstateL 
    FF_Lb2593taunu = FFEstateL 
    if sample == 'Data':
        weight = Weight_Data(tree)
    if sample == 'MISID':
        weight = Weight_MISID(tree)
    if sample == 'Combinatorial':
        weight = Weight_Combinatorial(tree,weighComb)
    if MCtype == 'MCfull':
        weight = getattr(tree,'Event_PIDCalibEffWeight')*getattr(tree,'w_LbCorr')
        if sample == 'Lb_Lcmunu' and FF_Lbmunu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
        if sample == 'Lb_Lctaunu' and FF_Lbtaunu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
        if sample == 'Lb_Lc2625munu' and FF_Lb2625munu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
        if sample == 'Lb_Lc2625taunu' and FF_Lb2625taunu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
        if sample == 'Lb_Lc2593munu' and FF_Lb2593munu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
        if sample == 'Lb_Lc2593taunu' and FF_Lb2593taunu==True:
            weight = weight*getattr(tree,'Event_FFcorr')
    return weight

def GetTemplateData(cut):
    h = r.TH3F('h_data_'+str(cut),"data",4,-2,14,10,0,2600,10,-2,14)
    for polarity in polarities:
        fname = folderD+'Lb_Data_'+polarity+'_preselected_full_sw.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_ISOLATION_BDT>cut:
                continue
            weight = GetWeight('Data',t,'')
            if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
    h.SetDirectory(0)
    return h

def GetTemplateMISID(cut):
    h = r.TH3F('h_MISID',"MISID",4,-2,14,10,0,2600,10,-2,14)
    for polarity in polarities:
        fname = folderMISID+'K_sample_'+polarity+'_full_sw_withCF.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_ISOLATION_BDT>cut:
                continue
            weight = GetWeight('MISID',t,'')
            if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6
,weight)
        fname = folderMISID+'Pi_sample_'+polarity+'_full_sw_withCF.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_ISOLATION_BDT>cut:
                continue
            weight = GetWeight('MISID',t,'')
            if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
    h.SetDirectory(0)
    return h
        
def GetTemplateCombinatorial(cut):
    h = r.TH3F('h_Combinatorial',"Combinatorial",4,-2,14,10,0,2600,10,-2,14)
    for polarity in polarities:
        fname = folderComb+'CombinatorialBkg_'+polarity+'_full.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            if t.Lb_ISOLATION_BDT>cut:
                continue
            weight = GetWeight('Combinatorial',t,'')
            if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
    h.SetDirectory(0)
    return h

def GetTemplateMCsample(mcsample, cut):
    if mcsample not in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        h = r.TH3F('h_'+mcsample,"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h.SetDirectory(0)
        for polarity in polarities:
            fname = folderMC+mcsample+'_'+polarity+'_full.root'
            fpresel_name = folderMC+mcsample+'_'+polarity+'_full_preselectionVars.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            fpresel = r.TFile(fpresel_name,'READ')
            tpresel = fpresel.Get('DecayTree')
            t.AddFriend(tpresel)
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if t.FinalSel!=1:
                    continue
                if t.Lb_ISOLATION_BDT>cut:
                    continue
                weight = GetWeight(mcsample,t,'MCfull')
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
        return h
    if mcsample in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        h_nominal = r.TH3F('h_'+mcsample,"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_2body = r.TH3F('h_'+mcsample+'-2body',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_mbody = r.TH3F('h_'+mcsample+'-mbody',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_mbody_1pl = r.TH3F('h_'+mcsample+'-mbody_1pl',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_mbody_1ml = r.TH3F('h_'+mcsample+'-mbody_1ml',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_mbody_1pq = r.TH3F('h_'+mcsample+'-mbody_1pq',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        h_mbody_1mq = r.TH3F('h_'+mcsample+'-mbody_1mq',"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)

        h_nominal.SetDirectory(0)
        h_2body.SetDirectory(0)
        h_mbody.SetDirectory(0)
        h_mbody_1pl.SetDirectory(0)
        h_mbody_1ml.SetDirectory(0)
        h_mbody_1pq.SetDirectory(0)
        h_mbody_1mq.SetDirectory(0)
        for polarity in polarities:
            fname = folderMC+'Lb_'+mcsample+'_'+polarity+'_full.root'
            fpresel_name = folderMC+'Lb_'+mcsample+'_'+polarity+'_full_preselectionVars.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            fpresel = r.TFile(fpresel_name,'READ')
            tpresel = fpresel.Get('DecayTree')
            t.AddFriend(tpresel)
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if t.FinalSel!=1:
                    continue
                if t.Lb_ISOLATION_BDT>cut:
                    continue
                weight = GetWeight(mcsample,t,'MCfull')
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    #Fill nominal histo (not used)
                    h_nominal.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    #Fill 2body histo
                    if t.twobody==True and t.mbody==False:
                        h_2body.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    #Fill 3body histo
                    if t.twobody==False and t.mbody==True:
                        h_mbody.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                        if t.w_mbody>-1000 and t.w_mbody!=1:
                            h_mbody_1pl.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_linear(t.w_mbody,1))
                            h_mbody_1ml.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_linear(t.w_mbody,-1))
                            h_mbody_1pq.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_quadratic(t.w_mbody,1))
                            h_mbody_1mq.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_quadratic(t.w_mbody,-1))

        h_mbody_1pl = ScaleHisto(h_mbody_1pl,h_mbody.Integral())
        h_mbody_1ml = ScaleHisto(h_mbody_1ml,h_mbody.Integral())
        h_mbody_1pq = ScaleHisto(h_mbody_1pq,h_mbody.Integral())
        h_mbody_1mq = ScaleHisto(h_mbody_1mq,h_mbody.Integral())
        return h_nominal

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

def PlotTemplateComparison(name,hold,hnew,norm):
    if norm:
        hold = ScaleHisto(hold,1)
        hnew = ScaleHisto(hnew,1)
    c = r.TCanvas('c_'+name,name,1500,500)
    c.Divide(3,1)
    c.cd(1)
    holdx = hold.ProjectionX()
    hnewx = hnew.ProjectionX()
    holdx.SetLineColor(r.kRed+2)
    hnewx.SetLineColor(r.kBlue+4)
    holdx.Draw('hist')
    hnewx.Draw('hist same')
    c.cd(2)
    holdy = hold.ProjectionY()
    hnewy = hnew.ProjectionY()
    holdy.SetLineColor(r.kRed+2)
    hnewy.SetLineColor(r.kBlue+4)
    holdy.Draw('hist')
    hnewy.Draw('hist same')
    c.cd(3)
    holdz = hold.ProjectionZ()
    hnewz = hnew.ProjectionZ()
    holdz.SetLineColor(r.kRed+2)
    hnewz.SetLineColor(r.kBlue+4)
    holdz.Draw('hist')
    hnewz.Draw('hist same')
    return c

h_Data_std = GetTemplateData(0.35)
h_Data_new = GetTemplateData(0.3)
c_data = PlotTemplateComparison('Data',h_Data_std,h_Data_new,0)
c_data.Draw()

h_MISID_std = GetTemplateMISID(0.35)
h_MISID_new = GetTemplateMISID(0.3)
c_MISID = PlotTemplateComparison('MISID',h_MISID_std,h_MISID_new,0)
c_MISID.Draw()

h_Combinatorial_std = GetTemplateCombinatorial(0.35)
h_Combinatorial_new = GetTemplateCombinatorial(0.3)
c_Combinatorial = PlotTemplateComparison('Combinatorial',h_Combinatorial_std,h_Combinatorial_new,0)
c_Combinatorial.Draw()

h_LbLcmunu_std = GetTemplateMCsample('Lb_Lcmunu',0.35)
h_LbLcmunu_new = GetTemplateMCsample('Lb_Lcmunu',0.3)
c_LbLcmunu = PlotTemplateComparison('LbLcmunu',h_LbLcmunu_std,h_LbLcmunu_new,0)
c_LbLcmunu.Draw()

h_LbLctaunu_std = GetTemplateMCsample('Lb_Lctaunu',0.35)
h_LbLctaunu_new = GetTemplateMCsample('Lb_Lctaunu',0.3)
c_LbLctaunu = PlotTemplateComparison('LbLctaunu',h_LbLctaunu_std,h_LbLctaunu_new,0)
c_LbLctaunu.Draw()

h_LbLcDs_std = GetTemplateMCsample('Lb_LcDs',0.35)
h_LbLcDs_new = GetTemplateMCsample('Lb_LcDs',0.3)
c_LbLcDs = PlotTemplateComparison('LbLcDs',h_LbLcDs_std,h_LbLcDs_new,0)
c_LbLcDs.Draw()

h_LbLc2625munu_std = GetTemplateMCsample('Lb_Lc2625munu',0.35)
h_LbLc2625munu_new = GetTemplateMCsample('Lb_Lc2625munu',0.3)
c_LbLc2625munu = PlotTemplateComparison('LbLc2625munu',h_LbLc2625munu_std,h_LbLc2625munu_new,0)
c_LbLc2625munu.Draw()

h_LbLc2593munu_std = GetTemplateMCsample('Lb_Lc2593munu',0.35)
h_LbLc2593munu_new = GetTemplateMCsample('Lb_Lc2593munu',0.3)
c_LbLc2593munu = PlotTemplateComparison('LbLc2593munu',h_LbLc2593munu_std,h_LbLc2593munu_new,0)
c_LbLc2593munu.Draw()

h_LbLc2625taunu_std = GetTemplateMCsample('Lb_Lc2625taunu',0.35)
h_LbLc2625taunu_new = GetTemplateMCsample('Lb_Lc2625taunu',0.3)
c_LbLc2625taunu = PlotTemplateComparison('LbLc2625taunu',h_LbLc2625taunu_std,h_LbLc2625taunu_new,0)
c_LbLc2625taunu.Draw()

h_LbLc2593taunu_std = GetTemplateMCsample('Lb_Lc2593taunu',0.35)
h_LbLc2593taunu_new = GetTemplateMCsample('Lb_Lc2593taunu',0.3)
c_LbLc2593taunu = PlotTemplateComparison('LbLc2593taunu',h_LbLc2593taunu_std,h_LbLc2593taunu_new,0)
c_LbLc2593taunu.Draw()

h_LbLc2625Ds_std = GetTemplateMCsample('Lb_Lc2625Ds',0.35)
h_LbLc2625Ds_new = GetTemplateMCsample('Lb_Lc2625Ds',0.3)
c_LbLc2625Ds = PlotTemplateComparison('LbLc2625Ds',h_LbLc2625Ds_std,h_LbLc2625Ds_new,0)
c_LbLc2625Ds.Draw()

h_LbLc2593Ds_std = GetTemplateMCsample('Lb_Lc2593Ds',0.35)
h_LbLc2593Ds_new = GetTemplateMCsample('Lb_Lc2593Ds',0.3)
c_LbLc2593Ds = PlotTemplateComparison('LbLc2593Ds',h_LbLc2593Ds_std,h_LbLc2593Ds_new,0)
c_LbLc2593Ds.Draw()

