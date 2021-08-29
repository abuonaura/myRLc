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



'''
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
            fname = folder+sample+'_'+polarity+'_full_preselected_'+suffix[category]+'.root'
            print(fname)
        if MCtype=='MCTrackerOnly':
            fname = folder+sample+'_'+polarity+'_preselected_'+suffix[category]+'.root'
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
'''



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
