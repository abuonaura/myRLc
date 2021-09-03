import ROOT as r
import os, sys
import argparse
from array import array

#Variables to change to include/exclude some corrections
suffixl = []
weighComb = False 
FFGstate  = True
FFEstateL = False
FFEstateH = True

if FFGstate:
    suffixl.append('_FFGstate')
if FFEstateL:
    suffixl.append('_FFEstateL')        #Form factors excited state light
if FFEstateH:
    suffixl.append('_FFEstateH')        #Form factors excited state heavy
if not FFGstate and not FFEstateL and not FFEstateH:
    suffixl.append('_NoFFcorr')

def init():
    ap = argparse.ArgumentParser(description='Create Input Histos for both isolated/K-enriched category')
    ap.add_argument('-c','--category', type=str, dest='category', default=None)
    ap.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    ap.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    args = ap.parse_args()
    return args

def PrintBinContent(h,template):
    print('------------------------')
    print('- Template: ', template)
    print('- NbinsX: {}, NbinsY: {}, NbinsZ{}:'.format(h.GetNbinsX(), h.GetNbinsY(), h.GetNbinsZ()))
    for i in range(h.GetNbinsX()):
        for j in range(h.GetNbinsY()):
            for k in range(h.GetNbinsZ()):
                print(' BinX: {}; BinY: {}; BinZ: {}; Content: {}'.format(i+1,j+1,k+1,h.GetBinContent(i+1,j+1,k+1)))
    return

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

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

def GetMCSampleNames(MCtype):
    mcsamples = {'MCfull':['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds'],'MCTrackerOnly':['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2880munu','Lb_Lc2765munu','B_Lcpbarmunu']}
    return mcsamples[MCtype]

def GetDataSampleNames():
    datasamples = ['Data','MISID','Combinatorial']
    return datasamples

def GetFolder(datatype,MCtype,sample):
    ddir = '/disk/lhcb_data2/RLcMuonic2016/'
    if datatype=='DATA':
        folder = {'Data':ddir+'Data/','MISID':ddir+'MISID/OppositeSign/',
                'Combinatorial':ddir+'CombinatorialBkg/'}
        return folder[sample]
    elif datatype=='MC':
        folder = {'MCfull':ddir+'MC_full_new/','MCTrackerOnly':ddir+'MC_TrackerOnly/'}
        return folder[MCtype]


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

double Weight_MultibodyCharm_linear(double w_mbody,double alpha1)
{
    double weight = 1+2*alpha1*w_mbody;
    return weight;
}

double Weight_MultibodyCharm_quadratic(double w_mbody,double alpha2)
{
    double weight = (1-alpha2)+8*alpha2*w_mbody*w_mbody;
    return weight;
}

double GetSigmas(double w)

'''

r.gInterpreter.Declare(func_weight)

def PutTogetherPolarityHistos(h,mcsample):
    h_new = h['MagUp']
    h_new.Add(h['MagDown'])
    h_new.SetTitle(mcsample)
    h_new.SetDirectory(0)
    return h_new

def GetTemplateData(category, cut):
    h = {}
    folder = GetFolder('DATA','','Data')
    for polarity in polarities:
        fname = folder+'Lb_Data_'+polarity+'_preselected_'+suffix[category]+'_sw.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("weight","GetWeightData(sw_sig)")
        h[polarity] = df0.Histo3D(("h_Data_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
    return h

def GetTemplateMISID(category):
    h = {}
    folder = GetFolder('DATA','','MISID')
    for polarity in polarities:
        fname1 = folder+'K_sample_'+polarity+'_'+suffix[category]+'_sw_withCF.root'
        fname2 = folder+'Pi_sample_'+polarity+'_'+suffix[category]+'_sw_withCF.root'
        fname = folder+'MISID_sample_'+polarity+'_'+suffix[category]+'_sw_withCF.root'
        os.system('hadd -f %s %s %s' %(fname, fname1, fname2))
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("weight","GetWeightMISID(w_recomu_CF,sw_sig)")
        h[polarity] = df0.Histo3D(("h_MISID_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
        column_name_vector = r.std.vector('string')()
        column_name_vector.push_back("weight")
        print("Writing the MISID tree with weights ...")
        df0.Snapshot("DecayTree",'MISID_weights_'+polariy+'_'+category+'.root',column_name_vector)
        print("Tree written.")
        os.system('rm %s' %fname)
    return h
        
def GetTemplateCombinatorial(category,weighComb):
    h = {}
    folder = GetFolder('DATA','','Combinatorial')
    for polarity in polarities:
        fname = folder+'CombinatorialBkg_'+polarity+'_'+suffix[category]+'.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        df0 = r.RDataFrame(t)
        df0 = df0.Define("HighMassCorrection",str(weighComb))
        df0 = df0.Define("weight","GetWeightCombinatorial(w_MISID,sw_sig,w_comb,HighMassCorrection)")
        h[polarity] = df0.Histo3D(("h_Combinatorial_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
        h[polarity].SetDirectory(0)
        column_name_vector = r.std.vector('string')()
        column_name_vector.push_back("weight")
        print("Writing the Combinatorial tree with weights ...")
        df0.Snapshot("DecayTree",'Comb_weights_'+polariy+'_'+category+'.root',column_name_vector)
        print("Tree written.")
    return h
    
        

def GetTemplateMCsample(category, mcsample, MCtype, cut):
    print(mcsample)
    folder = GetFolder('MC',MCtype,'')
    if mcsample not in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        h = {}
        h1 = {}
        for polarity in polarities:
            h[polarity] = r.TH3D()
            fname = folderMC+mcsample+'_'+polarity+'_'+suffix[category]+'.root'
            f = r.TFile(fname,'READ')
            t = f.Get('tupleout/DecayTree')
            df0 = r.RDataFrame(t)
            if mcsample=='Lb_Lcmunu' or mcsample=='Lb_Lctaunu':
                df0 = df0.Define("FFcorr",str(FFGstate))
            if mcsample=='Lb_Lc2593munu' or mcsample=='Lb_Lc2593taunu':
                df0 = df0.Define("FFcorr",str(FFEstateL))
            if mcsample=='Lb_Lc2625munu' or mcsample=='Lb_Lc2625taunu':
                df0 = df0.Define("FFcorr",str(FFEstateH))
            if mcsample in ["Lb_Lcmunu","Lb_Lctaunu","Lb_Lc2625munu","Lb_Lc2625taunu","Lb_Lc2593munu","Lb_Lc2593taunu"]:
                df0 = df0.Define("weight","GetMCweightWithFF(Event_PIDCalibEffWeight,w_LbCorr, Event_FFcorr,FFcorr)")
            else:
                df0 = df0.Define("weight","GetMCweightNoFF(Event_PIDCalibEffWeight,w_LbCorr")
            h[polarity] = df0.Histo3D(("h_"+mcsample+"_"+polarity+"_"+str(cut), "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            h[polarity].SetDirectory(0)
        #h1 = PutTogetherPolarityHistos(h,mcsample)
        return h
    if mcsample in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
        h_nominal, h_2body, h_mbody, h_mbody_1pl, h_mbody_1ml, h_mbody_1pq, h_mbody_1mq = {},{},{},{},{},{},{}
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
            df0 = df0.Define("weight","GetMCweightNoFF(Event_PIDCalibEffWeight,w_LbCorr)")
            h_nominal[polarity] = df0.Histo3D(("h_"+mcsample+"_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            df1 = df0.Filter("twobody==true && mbody==false")
            h_2body[polarity] = df1.Histo3D(("h_"+mcsample+"-2body_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            ##ADDITION of w_mbody>-1000 && w_mbody!=1 on 30.08.21
            df2 = df0.Filter("twobody==false && mbody==true && w_mbody>-1000 && w_mbody!=1")
            h_mbody[polarity] = df2.Histo3D(("h_"+mcsample+"-mbody_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight")
            df2.Define("weight_1ml","Weight_MultibodyCharm_linear(w_mbody,"+str(-1)+")")
            df2.Define("weight_1pl","Weight_MultibodyCharm_linear(w_mbody,"+str(1)+")")
            df2.Define("weight_1mq","Weight_MultibodyCharm_quadratic(w_mbody,"+str(-1)+")")
            df2.Define("weight_1pq","Weight_MultibodyCharm_quadratic(w_mbody,"+str(1)+")")
            h_mbody_1ml[polarity] = df2.Histo3D(("h_"+mcsample+"-mbody_1ml_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight*weight_1ml")
            h_mbody_1pl[polarity] = df2.Histo3D(("h_"+mcsample+"-mbody_1pl_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight*weight_1pl")
            h_mbody_1mq[polarity] = df2.Histo3D(("h_"+mcsample+"-mbody_1mq_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight*weight_1mq")
            h_mbody_1pq[polarity] = df2.Histo3D(("h_"+mcsample+"-mbody_1pq_"+polarity+"_"+category, "", 10, 0., 2600.,4, -2.E6, 14.E6,10, -2.E6, 14.E6),"FitVar_El_mLc","FitVar_q2_mLc","FitVar_Mmiss2_mLc","weight*weight_1pq")

            h_nominal[polarity].SetDirectory(0)
            h_2body[polarity].SetDirectory(0)
            h_mbody[polarity].SetDirectory(0)
            h_mbody_1pl[polarity].SetDirectory(0)
            h_mbody_1ml[polarity].SetDirectory(0)
            h_mbody_1pq[polarity].SetDirectory(0)
            h_mbody_1mq[polarity].SetDirectory(0)
        return h_nominal, h_2body, h_mbody, h_mbody_1pl, h_mbody_1ml, h_mbody_1pq, h_mbody_1mq


def GetIntegralAndSigmas(category,hMISID, hComb):
    integral_MISID = array('d',[0])
    integral_Comb  = array('d',[0])
    sigma_MISID    = array('d',[0])
    sigma_Comb     = array('d',[0])
    for polarity in polarities:
        integral_MISID[0] += hMISID[polarity].Integral()
        integral_Comb[0] += hComb[polarity].Integral()
        #extract the sigmas
        fMISID = r.TFile('MISID_weights_'+polariy+'_'+category+'.root','READ')
        tMISID = fMISID.Get('DecayTree')
        for i in range(tMISID.GetEntries()):
            tMISID.GetEntry(i)
            sigma_MISID[0]+=(tMISID.weight*tMISID.weight)
        fComb = r.TFile('Comb_weights_'+polariy+'_'+category+'.root','READ')
        tComb = fComb.Get('DecayTree')
        for i in range(tComb.GetEntries()):
            tComb.GetEntry(i)
            sigma_Comb[0]+=(tComb.weight*tComb.weight)
    sigma_MISID[0]=r.TMath.Sqrt(sigma_MISID[0])
    sigma_Comb[0] = r.TMath.Sqrt(sigma_Comb[0])
    return integral_MISID, integral_Comb, sigma_MISID, sigma_Comb
    
if __name__== "__main__":
    args = init()
    category = args.category
    if args.MCfull==True:
        MCtype = 'MCfull'
    if args.MCTO==True:
        MCtype = 'MCTrackerOnly'

    print('Category: ', category)
    print('MCtype: ', MCtype)

    if category==None:
        print('Please choose a category: Isolated or Kenriched')
        sys.exit('Aborting code')
    if MCtype==None:
        print('Please select MC type: Full or TrackerOnly')
        sys.exit('Aborting code')

    ofname = 'TemplateFiles/Histos_'+category+'_'+MCtype+'_LbCorr'
    for suffix in suffixl:
        print(suffix)
        ofname+=suffix
        print(ofname)
    ofname+='.root'
    print(ofname)

    outputFile = r.TFile.Open(ofname,"RECREATE")
    outputFile.SetCompressionAlgorithm(1)
    integral_MISID = array('d',[0])
    integral_Comb  = array('d',[0])
    sigma_MISID    = array('d',[0])
    sigma_Comb     = array('d',[0])
    tree = r.TTree('tree','tree')
    tree.Branch('integral_MISID',integral_MISID,'integral_MISID/D')
    tree.Branch('sigma_MISID',sigma_MISID,'sigma_MISID/D')
    tree.Branch('integral_Comb',integral_Comb,'integral_Comb/D')
    tree.Branch('sigma_Comb',sigma_Comb,'sigma_Comb/D')

    histos, sigma_MISID[0], sigma_Comb[0] = FillHistograms(category, MCtype)
    for h in histos:
        if h.GetName()=='h_MISID':
            integral_MISID[0] = h.Integral()
        if h.GetName()=='h_Combinatorial':
            integral_Comb[0] = h.Integral()
        h.SetDirectory(outputFile)
        SetNonNullBinContent(h)
    
    print('Expected yield MISID: '+str(integral_MISID)+'   sigma: '+str(sigma_MISID)) 
    print('Expected yield Comb: '+str(integral_Comb)+'   sigma: '+str(sigma_Comb)) 
    
    tree.Fill()

    outputFile.Write()
    outputFile.Close()
