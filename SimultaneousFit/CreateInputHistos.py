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
        folder = {'MCfull':ddir+'MC_full_new/','MCTrackerOnly':ddir+'MC_TrackerOnly'}
        return folder[MCtype]

def GetFileName(category,datatype,MCtype,sample,polarity):
    folder = GetFolder(datatype,MCtype,sample)
    suffix = {'Isolated':'iso','Kenriched':'Kenr','Lcpipi':'Lcpipi'}
    if datatype=='DATA':
        if sample=='Data':
            fname = folder+'Lb_'+sample+'_'+polarity+'_preselected_'+suffix[category]+'_sw.root'
        if sample=='MISID':
            fname = []
            fname.append(folder+'K_sample_'+polarity+'_'+suffix[category]+'_sw_withCF.root')
            fname.append(folder+'Pi_sample_'+polarity+'_'+suffix[category]+'_sw_withCF.root')
        if sample=='Combinatorial':
            fname = folder+'CombinatorialBkg_'+polarity+'_'+suffix[category]+'.root'
    if datatype=='MC':
        if MCtype=='MCfull':
            fname = folder+sample+'_'+polarity+'_full_preselected_'+suffix[category]+'_LbCorr.root'
            print(fname)
        if MCtype=='MCTrackerOnly':
            fname = folder+sample+'_'+polarity+'_preselected_'+suffix[category]+'_LbCorr.root'
    return fname

def GetFileListPerSample(category,datatype,MCtype,sample):
    polarities = ['MagUp','MagDown']
    flist = [] 
    for polarity in polarities:
        fname = GetFileName(category,datatype,MCtype,sample,polarity)
        print (fname)
        if sample!='MISID':
            flist.append(fname)
        else:
            flist.extend(fname)
    print(flist)
    return flist


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
    if sample in GetMCSampleNames(MCtype):
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


def FillHistograms(category,MCtype):
    histos = []
    sigma_MISID = array('d',[0])
    sigma_Comb = array('d',[0])
    #HISTOGRAMS FOR DATA BASED TEMPLATES
    datasamples = GetDataSampleNames()
    for datasample in datasamples:
        print(datasample)
        filelist = GetFileListPerSample(category,'DATA','',datasample)
        print(filelist)
        h = r.TH3F('h_'+datasample,"qem_"+datasample,4,-2,14,10,0,2600,10,-2,14)
        h.SetDirectory(0)
        for fname in filelist:
            f = r.TFile(fname,'READ')
            t = f.Get('DecayTree')
            for i in range(t.GetEntries()):
                t.GetEntry(i)
                weight = GetWeight(datasample,t,MCtype)
                if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                    h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    if datasample=='MISID':
                        sigma_MISID[0]+=weight*weight
                    if datasample=='Combinatorial':
                        sigma_Comb[0]+=weight*weight
        histos.append(h)
    sigma_MISID[0]=r.TMath.Sqrt(sigma_MISID[0])
    sigma_Comb[0] = r.TMath.Sqrt(sigma_Comb[0])

    #HISTOGRAMS FOR MC BASED TEMPLATES
    mcsamples = GetMCSampleNames(MCtype)
    print(mcsamples)
    for mcsample in mcsamples:
        print(mcsample)
        filelist = GetFileListPerSample(category,'MC',MCtype,mcsample)
        print(filelist)
        if mcsample not in ['Lb_LcDs','Lb_Lc2625Ds','Lb_Lc2593Ds']:
            h = r.TH3F('h_'+mcsample,"qem_"+mcsample,4,-2,14,10,0,2600,10,-2,14)
            h.SetDirectory(0)
            for fname in filelist:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight = GetWeight(mcsample,t,MCtype)
                    if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                        h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
            histos.append(h)
            print(histos)
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
            for fname in filelist:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight = GetWeight(mcsample,t,MCtype)
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
            histos.append(h_nominal)
            histos.append(h_2body)
            histos.append(h_mbody)
            histos.append(h_mbody_1pl)
            histos.append(h_mbody_1ml)
            histos.append(h_mbody_1pq)
            histos.append(h_mbody_1mq)
    return histos, sigma_MISID[0], sigma_Comb[0]


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
