'''
Author: Annarita Buonaura
Date: November 2019

Description: Script to produce templates for fitting

How to run:
    - python CreateInputHistos.py -c Isolated (Kenriched) --MCfull (--MCTrackerOnly)
'''

import ROOT as r
import os, sys
import argparse

def init():
    ap = argparse.ArgumentParser(description='Create Input Histos for both isolated/K-enriched category')
    ap.add_argument('-c','--category', type=str, dest='category', default=None)
    ap.add_argument('--MCfull',dest='MCfull', help="Process MC full simulation samples", required=False, default=False, action='store_true')
    ap.add_argument('--MCTrackerOnly',dest='MCTO', help="Process MC TrackerOnly simulation samples", required=False, default=False, action='store_true')
    args = ap.parse_args()
    return args

def GetFileList(category,MCtype):
    if category=='Isolated':
        if MCtype=='MCfull':
            startfiles={'data':['Data/Lb_Data_MagUp_preselected_iso_sw.root','Data/Lb_Data_MagDown_preselected_iso_sw.root'],
                        'MISID':['MISID/OppositeSign/K_sample_MagDown_iso_sw_withCF.root','MISID/OppositeSign/K_sample_MagUp_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_iso_sw_withCF.root'],
                        'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_iso.root','CombinatorialBkg/CombinatorialBkg_MagDown_iso.root'],
                        'mu':['MC_full_new/Lb_Lcmunu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lcmunu_MagDown_full_preselected_iso_LbCorr.root'],
                        'tau':['MC_full_new/Lb_Lctaunu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lctaunu_MagDown_full_preselected_iso_LbCorr.root'],
                        '2charm':['MC_full_new/Lb_LcDs_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_LcDs_MagDown_full_preselected_iso_LbCorr.root'],
                        'starDs':['MC_full_new/Lb_Lc2625Ds_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2625Ds_MagDown_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagDown_full_preselected_iso_LbCorr.root'],
                        'starmu':['MC_full_new/Lb_Lc2625munu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2625munu_MagDown_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagDown_full_preselected_iso_LbCorr.root'],
                        'startau':['MC_full_new/Lb_Lc2625taunu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2625taunu_MagDown_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagUp_full_preselected_iso_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagDown_full_preselected_iso_LbCorr.root']}
        if MCtype=='MCTrackerOnly':
            startfiles={'data':['Data/Lb_Data_MagUp_preselected_iso_sw.root','Data/Lb_Data_MagDown_preselected_iso_sw.root'],
                        'MISID':['MISID/OppositeSign/K_sample_MagDown_iso_sw_withCF.root','MISID/OppositeSign/K_sample_MagUp_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_iso_sw_withCF.root'],
                        'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_iso.root','CombinatorialBkg/CombinatorialBkg_MagDown_iso.root'],
                        'mu':['MC_TrackerOnly/Lb_Lcmunu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lcmunu_MagDown_preselected_iso_LbCorr.root'],
                        'tau':['MC_TrackerOnly/Lb_Lctaunu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lctaunu_MagDown_preselected_iso_LbCorr.root'],
                        '2charm':['MC_TrackerOnly/Lb_LcDs_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_LcDs_MagDown_preselected_iso_LbCorr.root'],
                        'starDs':['MC_TrackerOnly/Lb_Lc2625Ds_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2625Ds_MagDown_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593Ds_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593Ds_MagDown_preselected_iso_LbCorr.root'],
                        'starmu':['MC_TrackerOnly/Lb_Lc2625munu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2625munu_MagDown_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593munu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593munu_MagDown_preselected_iso_LbCorr.root'],
                        'startau':['MC_TrackerOnly/Lb_Lc2625taunu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2625taunu_MagDown_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593taunu_MagUp_preselected_iso_LbCorr.root','MC_TrackerOnly/Lb_Lc2593taunu_MagDown_preselected_iso_LbCorr.root']}

    
    if category=='Kenriched':
        if MCtype=='MCfull':
            startfiles ={'data':['Data/Lb_Data_MagUp_preselected_Kenr_sw.root','Data/Lb_Data_MagDown_preselected_Kenr_sw.root'],
            'MISID':['MISID/OppositeSign/K_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/K_sample_MagDown_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_Kenr_sw_withCF.root'],
            'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_Kenr.root','CombinatorialBkg/CombinatorialBkg_MagDown_Kenr.root'],
            'mu':['MC_full_new/Lb_Lcmunu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lcmunu_MagDown_full_preselected_Kenr_LbCorr.root'],
            'tau':['MC_full_new/Lb_Lctaunu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lctaunu_MagDown_full_preselected_Kenr_LbCorr.root'],
            '2charm':['MC_full_new/Lb_LcDs_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_LcDs_MagDown_full_preselected_Kenr_LbCorr.root'],
            'starmu':['MC_full_new/Lb_Lc2625munu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2625munu_MagDown_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagDown_full_preselected_Kenr_LbCorr.root'],
            'starDs':['MC_full_new/Lb_Lc2625Ds_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2625Ds_MagDown_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagDown_full_preselected_Kenr_LbCorr.root'],
            'startau':['MC_full_new/Lb_Lc2625taunu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagUp_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2625taunu_MagDown_full_preselected_Kenr_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagDown_full_preselected_Kenr_LbCorr.root']
            }
        if MCtype=='MCTrackerOnly':
            startfiles ={'data':['Data/Lb_Data_MagUp_preselected_Kenr_sw.root','Data/Lb_Data_MagDown_preselected_Kenr_sw.root'],
            'MISID':['MISID/OppositeSign/K_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/K_sample_MagDown_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_Kenr_sw_withCF.root'],
            'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_Kenr.root','CombinatorialBkg/CombinatorialBkg_MagDown_Kenr.root'],
            'mu':['MC_TrackerOnly/Lb_Lcmunu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lcmunu_MagDown_preselected_Kenr_LbCorr.root'],
            'tau':['MC_TrackerOnly/Lb_Lctaunu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lctaunu_MagDown_preselected_Kenr_LbCorr.root'],
            '2charm':['MC_TrackerOnly/Lb_LcDs_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_LcDs_MagDown_preselected_Kenr_LbCorr.root'],
            'starmu':['MC_TrackerOnly/Lb_Lc2625munu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593munu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2625munu_MagDown_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593munu_MagDown_preselected_Kenr_LbCorr.root'],
            'starDs':['MC_TrackerOnly/Lb_Lc2625Ds_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593Ds_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2625Ds_MagDown_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593Ds_MagDown_preselected_Kenr_LbCorr.root'],
            'startau':['MC_TrackerOnly/Lb_Lc2625taunu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593taunu_MagUp_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2625taunu_MagDown_preselected_Kenr_LbCorr.root','MC_TrackerOnly/Lb_Lc2593taunu_MagDown_preselected_Kenr_LbCorr.root']
            }
    if category=='Lcpipi':
        if MCtype=='MCfull':
            startfiles ={'data':['Data/Lb_Data_MagUp_preselected_Lcpipi_sw.root','Data/Lb_Data_MagDown_preselected_Lcpipi_sw.root'],
            'MISID':['MISID/OppositeSign/K_sample_MagUp_Lcpipi_sw_withCF.root','MISID/OppositeSign/K_sample_MagDown_Lcpipi_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_Lcpipi_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_Lcpipi_sw_withCF.root'],
            'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_Lcpipi.root','CombinatorialBkg/CombinatorialBkg_MagDown_Lcpipi.root'],
            'mu':['MC_full_new/Lb_Lcmunu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lcmunu_MagDown_full_preselected_Lcpipi_LbCorr.root'],
            'tau':['MC_full_new/Lb_Lctaunu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lctaunu_MagDown_full_preselected_Lcpipi_LbCorr.root'],
            '2charm':['MC_full_new/Lb_LcDs_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_LcDs_MagDown_full_preselected_Lcpipi_LbCorr.root'],
            'starmu':['MC_full_new/Lb_Lc2625munu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2625munu_MagDown_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593munu_MagDown_full_preselected_Lcpipi_LbCorr.root'],
            'starDs':['MC_full_new/Lb_Lc2625Ds_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2625Ds_MagDown_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593Ds_MagDown_full_preselected_Lcpipi_LbCorr.root'],
            'startau':['MC_full_new/Lb_Lc2625taunu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagUp_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2625taunu_MagDown_full_preselected_Lcpipi_LbCorr.root','MC_full_new/Lb_Lc2593taunu_MagDown_full_preselected_Lcpipi_LbCorr.root']
            }

    return startfiles

def GetStartRootFile(startfiles,sample):
    datadir = '/disk/lhcb_data2/RLcMuonic2016/'
    rootfiles = []
    for iF in startfiles[sample]:
        iF = datadir+iF
        rootfiles.append(iF)
    return rootfiles

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

def ApplyIsolationMC(tree):
    if getattr(tree,'Lb_ISOLATION_BDT')<0.35:
        flag = True
    else:
        flag = False
    return flag

def TruthMatch(tree):
    if getattr(tree,'Lc_BKGCAT')<30 and getattr(tree,'Lb_BKGCAT')<50:
        flag = True
    else:
        flag = False
    return flag

def Weight_MultibodyCharm_linear(tree,alpha1):
    weight = 1+2*alpha1*getattr(tree,'w_mbody')
    return weight

def Weight_MultibodyCharm_quadratic(tree,alpha2):
    weight = (1-alpha2)+8*alpha2*getattr(tree,'w_mbody')*getattr(tree,'w_mbody')
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

def FillHistograms(category, MCtype, nbinsq2, nbinsEl, nbinsMM):
    sigma_MISID = 0
    sigma_Comb = 0
    histos = []
    startfiles = GetFileList(category,MCtype)

    n2charm_2body=0
    n2charm_mbody=0

    for sample in startfiles.keys():
        print(sample)
        files = GetStartRootFile(startfiles,sample)
        print(files)
        if sample!='2charm' and sample!='starDs':
            h = r.TH3F('h_'+sample,"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h.SetDirectory(0)
            for fname in files:
                integral = 0
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight=0
                    if sample == 'data':
                        weight = Weight_Data(t)
                    if sample == 'MISID':
                        if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                            weight = Weight_MISID(t)
                            sigma_MISID = sigma_MISID + weight*weight
                    if sample == 'Combinatorial':
                        weight = Weight_Combinatorial(t,False)
                        if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                            integral+=weight
                            sigma_Comb = sigma_Comb+ weight*weight
                    if sample == 'mu' or sample=='tau' or sample=='starmu' or sample=='startau':
                        weight=t.Event_PIDCalibEffWeight*t.w_LbCorr
                    h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
            histos.append(h)
            if sample=='MISID' or sample=='Combinatorial':
                print(' Values for eventual Gaussian constraint: ')
            if sample=='MISID':
                sigma_MISID = r.TMath.Sqrt(sigma_MISID)
                print('Sample: ', sample)
                print ('   Ex. yield: ', h.Integral())
                print ('   sigma: ', sigma_MISID)
            if sample=='Combinatorial':
                sigma_Comb = r.TMath.Sqrt(sigma_Comb)
                print('Sample: ', sample)
                print ('   Ex. yield: ', h.Integral())
                print ('   sigma: ', sigma_Comb)

        if sample=='2charm' or sample=='starDs':
            h_nominal = r.TH3F('h_'+sample,"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_2body = r.TH3F('h_'+sample+'-2body',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_mbody = r.TH3F('h_'+sample+'-mbody',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_mbody_1pl = r.TH3F('h_'+sample+'-mbody_1pl',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_mbody_1ml = r.TH3F('h_'+sample+'-mbody_1ml',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_mbody_1pq = r.TH3F('h_'+sample+'-mbody_1pq',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)
            h_mbody_1mq = r.TH3F('h_'+sample+'-mbody_1mq',"qem_"+sample,nbinsq2,-2,14,nbinsEl,0,2600,nbinsMM,-2,14)

            h_nominal.SetDirectory(0)
            h_2body.SetDirectory(0)
            h_mbody.SetDirectory(0)
            h_mbody_1pl.SetDirectory(0)
            h_mbody_1ml.SetDirectory(0)
            h_mbody_1pq.SetDirectory(0)
            h_mbody_1mq.SetDirectory(0)

            for fname in files:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight = t.Event_PIDCalibEffWeight*t.w_LbCorr
                    h_nominal.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    if t.twobody==True and t.mbody==False:
                        n2charm_2body+=1
                        h_2body.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                    if t.twobody==False and t.mbody==True:
                        n2charm_mbody+=1
                        h_mbody.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
                        if t.w_mbody>-1000 and t.w_mbody!=1:
                            h_mbody_1pl.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_linear(t,1))
                            h_mbody_1ml.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_linear(t,-1))
                            h_mbody_1pq.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_quadratic(t,1))
                            h_mbody_1mq.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight*Weight_MultibodyCharm_quadratic(t,-1))
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

            print("n2charm_2body: ", n2charm_2body)
            print("n2charm_mbody: ", n2charm_mbody)
    return histos

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

    outfile_name = 'RootFiles/Histos_'+category+'_'+MCtype+'_morebins.root'
    outputFile = r.TFile.Open(outfile_name,"RECREATE")
    outputFile.SetCompressionAlgorithm(1)
    histos = FillHistograms(category, MCtype,10,20,20)
    for h in histos:
        h.SetDirectory(outputFile)
        SetNonNullBinContent(h)
    outputFile.Write()
    outputFile.Close()

