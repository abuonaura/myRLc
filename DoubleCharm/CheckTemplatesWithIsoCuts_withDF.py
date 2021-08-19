import ROOT as r
from ROOT import TFile, TTree, TH1F, TCanvas, TMath, TLegend
from array import array
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys, os, glob

folderD = '/disk/lhcb_data2/RLcMuonic2016/Data/'
folderMISID = '/disk/lhcb_data2/RLcMuonic2016/MISID/OppositeSign/'
folderComb = '/disk/lhcb_data2/RLcMuonic2016/CombinatorialBkg/'
folderMC = '/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/'

datasamples = ['Data','MISID','Combinatorial']
mcsamples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds']
#polarities = ['MagUp','MagDown']
polarities = ['MagUp']

#Variables to change to include/exclude some corrections
suffixl = []
weighComb = 0
FFGstate  = 0
FFEstateL = 0
FFEstateH = 0

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

def GetWeight(sample,tree):
    if sample == 'Data':
        weight = Weight_Data(tree)
    if sample == 'MISID':
        weight = Weight_MISID(tree)
    if sample == 'Combinatorial':
        weight = Weight_Combinatorial(tree,weighComb)
    return weight

func_weight = '''
double GetWeightData(double sw_sig)
{
    return sw_sig;
}

double GetWeightMISID(double w_recomu_CF,double sw_sig)
{
    double weight = w_recomu_CF*10*sw_sig;
    return weight;
}

double GetWeightCombinatorial(double w_MISID, double sw_sig, double w_comb, bool HighMassCorrection)
{
    double weight = w_MISID*sw_sig;
    if(HighMassCorrection)
        weight =  w_MISID*sw_sig*w_comb;
    return weight;
}

double GetMCweightWithFF(double Event_PIDCalibEffWeight, double w_LbCorr, double Event_FFcorr, bool FFcorr)
{
    double weight = Event_PIDCalibEffWeight*w_LbCorr;
    if(FFcorr)
        weight = weight*Event_FFcorr;
    return weight;
}
double GetMCweightNoFF(double Event_PIDCalibEffWeight, double w_LbCorr)
{
    double weight = Event_PIDCalibEffWeight*w_LbCorr;
    return weight;
}
'''

r.gInterpreter.Declare(func_weight)

def GetTemplateData(cut):
    h = {}
    for polarity in polarities:
        fname = folderD+'Lb_Data_'+polarity+'_preselected_full_sw.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("cut",str(cut))
        df1 = df0.Filter("Lb_ISOLATION_BDT<=cut")
        df1 = df1.Define("weight","GetWeightData(sw_sig)")
        h[polarity] = df1.Histo3D(("h_Data_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
    return h['MagUp']

def GetTemplateMISID(cut):
    h = {}
    for polarity in polarities:
        fname1 = folderMISID+'K_sample_'+polarity+'_full_sw_withCF.root'
        fname2 = folderMISID+'Pi_sample_'+polarity+'_full_sw_withCF.root'
        fname = folderMISID+'MISID_sample_'+polarity+'_full_sw_withCF.root'
        os.system('hadd -f %s %s %s' %(fname, fname1, fname2))
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("cut",str(cut))
        df1 = df0.Filter("Lb_ISOLATION_BDT<=cut")
        df1 = df1.Define("weight","GetWeightMISID(w_recomu_CF,sw_sig)")
        h[polarity] = df1.Histo3D(("h_MISID_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
        os.system('rm %s' %fname)
    return h['MagUp']
        
def GetTemplateCombinatorial(cut):
    h = {}
    for polarity in polarities:
        fname = folderComb+'CombinatorialBkg_'+polarity+'_full.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("cut",str(cut))
        df0 = df0.Define("HighMassCorrection",str(weighComb))
        df1 = df0.Filter("Lb_ISOLATION_BDT<=cut")
        df1 = df1.Define("weight","GetWeightCombinatorial(w_MISID,sw_sig,w_comb,HighMassCorrection)")
        h[polarity] = df1.Histo3D(("h_MISID_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
    return h['MagUp']

def PutTogetherPolarityHistos(h,mcsample):
    h_new = h['MagUp']
    h_new.Add(h['MagDown'])
    h_new.SetTitle(mcsample)
    h_new.SetDirectory(0)
    return h_new

def GetTemplateMCsample(mcsample, cut):
    print(mcsample)
    if mcsample not in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        #h = r.TH3F('h_'+mcsample,"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
        #h.SetDirectory(0)
        h = {}
        h1 = {}
        for polarity in polarities:
            h[polarity] = r.TH3D()
            fname = folderMC+mcsample+'_'+polarity+'_full.root'
            fpresel_name = folderMC+mcsample+'_'+polarity+'_full_preselectionVars.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            t.SetBranchStatus("*",0)
            t.SetBranchStatus('Lb_ISOLATION_BDT',1)
            t.SetBranchStatus('FitVar_q2_mLc',1)
            t.SetBranchStatus('FitVar_El_mLc',1)
            t.SetBranchStatus('FitVar_Mmiss2_mLc',1)
            fpresel = r.TFile(fpresel_name,'READ')
            tpresel = fpresel.Get('DecayTree')
            t.AddFriend(tpresel)
            df0 = r.RDataFrame(t)
            df0 = df0.Define("cut",str(cut))
            df1 = df0.Filter("FinalSel==1&&Lb_ISOLATION_BDT<=cut")
            #df1 = df1.Define("sample",str(mcsample))
            if mcsample=='Lb_Lcmunu' or mcsample=='Lb_Lctaunu':
                df1 = df1.Define("FFcorr",str(FFGstate))
            if mcsample=='Lb_Lc2593munu' or mcsample=='Lb_Lc2593taunu':
                df1 = df1.Define("FFcorr",str(FFEstateL))
            if mcsample=='Lb_Lc2625munu' or mcsample=='Lb_Lc2625taunu':
                df1 = df1.Define("FFcorr",str(FFEstateH))
            if mcsample in ["Lb_Lcmunu","Lb_Lctaunu","Lb_Lc2625munu","Lb_Lc2625taunu","Lb_Lc2593munu","Lb_Lc2593taunu"]:
                df1 = df1.Define("weight","GetMCweightWithFF(Event_PIDCalibEffWeight,w_LbCorr, Event_FFcorr,FFcorr)")
            else:
                df1 = df1.Define("weight","GetMCweightNoFF(Event_PIDCalibEffWeight,w_LbCorr")
            h[polarity] = df1.Histo3D(("h_"+mcsample+"_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            '''
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                if t.FinalSel!=1 and t.Lb_ISOLATION_BDT>cut:
                    continue
                if (i%100000==0):
                    print('Processed evt ',i)
                weight = GetWeight(mcsample,t,'MCfull')
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
            '''
            h[polarity].SetDirectory(0)
        #h1 = PutTogetherPolarityHistos(h,mcsample)
        return h['MagUp']
    if mcsample in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        h_nominal, h_2body, h_mbody = {},{},{}
        for polarity in polarities:
            fname = folderMC+mcsample+'_'+polarity+'_full.root'
            fpresel_name = folderMC+mcsample+'_'+polarity+'_full_preselectionVars.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            t.SetBranchStatus("*",0)
            t.SetBranchStatus('Lb_ISOLATION_BDT',1)
            t.SetBranchStatus('FitVar_q2_mLc',1)
            t.SetBranchStatus('FitVar_El_mLc',1)
            t.SetBranchStatus('FitVar_Mmiss2_mLc',1)
            fpresel = r.TFile(fpresel_name,'READ')
            tpresel = fpresel.Get('DecayTree')
            t.AddFriend(tpresel)
            df0 = r.RDataFrame(t)
            df0 = df0.Define("cut",str(cut))
            df1 = df0.Filter("FinalSel==1&&Lb_ISOLATION_BDT<=cut")
            df1 = df1.Define("weight","GetMCweightNoFF(Event_PIDCalibEffWeight,w_LbCorr)")
            h_nominal[polarity] = df1.Histo3D(("h_"+mcsample+"_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            df2 = df1.Filter("twobody==true && mbody==false")
            h_2body[polarity] = df2.Histo3D(("h_2body_"+mcsample+"_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            df3 = df1.Filter("twobody==false && mbody==true")
            h_mbody[polarity] = df3.Histo3D(("h_mbody_"+mcsample+"_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            h_nominal[polarity].SetDirectory(0)
            h_2body[polarity].SetDirectory(0)
            h_mbody[polarity].SetDirectory(0)
        return h_nominal['MagUp']

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
c_data = PlotTemplateComparison('Data',h_Data_std,h_Data_new,1)
c_data.Draw()
c_data.SaveAs("plots/ComparisonData.png")

h_MISID_std = GetTemplateMISID(0.35)
h_MISID_new = GetTemplateMISID(0.3)
c_MISID = PlotTemplateComparison('MISID',h_MISID_std,h_MISID_new,1)
c_MISID.Draw()
c_MISID.SaveAs("plots/ComparisonMISID.png")

h_Combinatorial_std = GetTemplateCombinatorial(0.35)
h_Combinatorial_new = GetTemplateCombinatorial(0.3)
c_Combinatorial = PlotTemplateComparison('Combinatorial',h_Combinatorial_std,h_Combinatorial_new,1)
c_Combinatorial.Draw()
c_Combinatorial.SaveAs("plots/ComparisonCombinatorial.png")

h_LbLcmunu_std = GetTemplateMCsample('Lb_Lcmunu',0.35)
h_LbLcmunu_new = GetTemplateMCsample('Lb_Lcmunu',0.3)
c_LbLcmunu = PlotTemplateComparison('LbLcmunu',h_LbLcmunu_std,h_LbLcmunu_new,1)
c_LbLcmunu.Draw()
c_LbLcmunu.SaveAs("plots/ComparisonLbLcmunu.png")

h_LbLctaunu_std = GetTemplateMCsample('Lb_Lctaunu',0.35)
h_LbLctaunu_new = GetTemplateMCsample('Lb_Lctaunu',0.3)
c_LbLctaunu = PlotTemplateComparison('LbLctaunu',h_LbLctaunu_std,h_LbLctaunu_new,1)
c_LbLctaunu.Draw()
c_LbLctaunu.SaveAs("plots/ComparisonLbLctaunu.png")

h_LbLcDs_std = GetTemplateMCsample('Lb_LcDs',0.35)
h_LbLcDs_new = GetTemplateMCsample('Lb_LcDs',0.3)
c_LbLcDs = PlotTemplateComparison('LbLcDs',h_LbLcDs_std,h_LbLcDs_new,1)
c_LbLcDs.Draw()
c_LbLcDs.SaveAs("plots/ComparisonLbLcDs.png")

h_LbLc2625munu_std = GetTemplateMCsample('Lb_Lc2625munu',0.35)
h_LbLc2625munu_new = GetTemplateMCsample('Lb_Lc2625munu',0.3)
c_LbLc2625munu = PlotTemplateComparison('LbLc2625munu',h_LbLc2625munu_std,h_LbLc2625munu_new,1)
c_LbLc2625munu.Draw()
c_LbLc2625munu.SaveAs("plots/ComparisonLbLc2625munu.png")

h_LbLc2593munu_std = GetTemplateMCsample('Lb_Lc2593munu',0.35)
h_LbLc2593munu_new = GetTemplateMCsample('Lb_Lc2593munu',0.3)
c_LbLc2593munu = PlotTemplateComparison('LbLc2593munu',h_LbLc2593munu_std,h_LbLc2593munu_new,1)
c_LbLc2593munu.Draw()
c_LbLc2593munu.SaveAs("plots/ComparisonLbLc2593munu.png")

h_LbLc2625taunu_std = GetTemplateMCsample('Lb_Lc2625taunu',0.35)
h_LbLc2625taunu_new = GetTemplateMCsample('Lb_Lc2625taunu',0.3)
c_LbLc2625taunu = PlotTemplateComparison('LbLc2625taunu',h_LbLc2625taunu_std,h_LbLc2625taunu_new,1)
c_LbLc2625taunu.Draw()
c_LbLc2625taunu.SaveAs("plots/ComparisonLbLc2625taunu.png")

h_LbLc2593taunu_std = GetTemplateMCsample('Lb_Lc2593taunu',0.35)
h_LbLc2593taunu_new = GetTemplateMCsample('Lb_Lc2593taunu',0.3)
c_LbLc2593taunu = PlotTemplateComparison('LbLc2593taunu',h_LbLc2593taunu_std,h_LbLc2593taunu_new,1)
c_LbLc2593taunu.Draw()
c_LbLc2593taunu.SaveAs("plots/ComparisonLbLc2593taunu.png")

h_LbLc2625Ds_std = GetTemplateMCsample('Lb_Lc2625Ds',0.35)
h_LbLc2625Ds_new = GetTemplateMCsample('Lb_Lc2625Ds',0.3)
c_LbLc2625Ds = PlotTemplateComparison('LbLc2625Ds',h_LbLc2625Ds_std,h_LbLc2625Ds_new,1)
c_LbLc2625Ds.Draw()
c_LbLc2625Ds.SaveAs("plots/ComparisonLbLc2625Ds.png")

h_LbLc2593Ds_std = GetTemplateMCsample('Lb_Lc2593Ds',0.35)
h_LbLc2593Ds_new = GetTemplateMCsample('Lb_Lc2593Ds',0.3)
c_LbLc2593Ds = PlotTemplateComparison('LbLc2593Ds',h_LbLc2593Ds_std,h_LbLc2593Ds_new,1)
c_LbLc2593Ds.Draw()
c_LbLc2593Ds.SaveAs("plots/ComparisonLbLc2593Ds.png")

fout = r.TFile('TemplatesForIsoCutsComparison.root','RECREATE')
h_Data_std.Write()
h_Data_new.Write()
h_MISID_std.Write()
h_MISID_new.Write()
h_Combinatorial_std.Write()
h_Combinatorial_new.Write()
h_LbLcmunu_std.Write()
h_LbLcmunu_new.Write()
h_LbLctaunu_std.Write()
h_LbLctaunu_new.Write()
h_LbLcDs_std.Write()
h_LbLcDs_new.Write()
h_LbLc2593munu_std.Write()
h_LbLc2593munu_new.Write()
h_LbLc2593taunu_std.Write()
h_LbLc2593taunu_new.Write()
h_LbLc2593Ds_std.Write()
h_LbLc2593Ds_new.Write()
h_LbLc2625munu_std.Write()
h_LbLc2625munu_new.Write()
h_LbLc2625taunu_std.Write()
h_LbLc2625taunu_new.Write()
h_LbLc2625Ds_std.Write()
h_LbLc2625Ds_new.Write()
fout.Close()








