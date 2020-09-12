import ROOT as r

def GetFileList(category):
    if category=='Kenriched':
        startfiles ={'data':['Data/Lb_Data_MagUp_preselected_Kenr_sw.root','Data/Lb_Data_MagDown_preselected_Kenr_sw.root'],
            'MISID':['MISID/OppositeSign/K_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/K_sample_MagDown_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_Kenr_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_Kenr_sw_withCF.root'],
            'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_Kenr.root','CombinatorialBkg/CombinatorialBkg_MagDown_Kenr.root'],
            'mu':['MC_full_new/Lb_Lcmunu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lcmunu_MagDown_full_preselected_Kenr.root'],
            'tau':['MC_full_new/Lb_Lctaunu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lctaunu_MagDown_full_preselected_Kenr.root'],
            '2charm':['MC_full_new/Lb_LcDs_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_LcDs_MagDown_full_preselected_Kenr.root'],
            'starmu':['MC_full_new/Lb_Lc2625munu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593munu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2625munu_MagDown_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593munu_MagDown_full_preselected_Kenr.root'],
            'starDs':['MC_full_new/Lb_Lc2625Ds_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593Ds_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2625Ds_MagDown_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593Ds_MagDown_full_preselected_Kenr.root'],
            'startau':['MC_full_new/Lb_Lc2625taunu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593taunu_MagUp_full_preselected_Kenr.root','MC_full_new/Lb_Lc2625taunu_MagDown_full_preselected_Kenr.root','MC_full_new/Lb_Lc2593taunu_MagDown_full_preselected_Kenr.root']
            }
    if category=='Isolated':
        startfiles={'data':['Data/Lb_Data_MagUp_preselected_iso_sw.root','Data/Lb_Data_MagDown_preselected_iso_sw.root'],
                        'MISID':['MISID/OppositeSign/K_sample_MagDown_iso_sw_withCF.root'
,'MISID/OppositeSign/K_sample_MagUp_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagDown_iso_sw_withCF.root','MISID/OppositeSign/Pi_sample_MagUp_iso_sw_withCF.root'],
                        'Combinatorial': ['CombinatorialBkg/CombinatorialBkg_MagUp_iso.root','CombinatorialBkg/CombinatorialBkg_MagDown_iso.root'],
                        'mu':['MC_full_new/Lb_Lcmunu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lcmunu_MagDown_full_preselected_iso.root'],
                        'tau':['MC_full_new/Lb_Lctaunu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lctaunu_MagDown_full_preselected_iso.root'],
                        '2charm':['MC_full_new/Lb_LcDs_MagUp_full_preselected_iso.root','MC_full_new/Lb_LcDs_MagDown_full_preselected_iso.root'],
                        'starDs':['MC_full_new/Lb_Lc2625Ds_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2625Ds_MagDown_full_preselected_iso.root','MC_full_new/Lb_Lc2593Ds_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2593Ds_MagDown_full_preselected_iso.root'],
                        'starmu':['MC_full_new/Lb_Lc2625munu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2625munu_MagDown_full_preselected_iso.root','MC_full_new/Lb_Lc2593munu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2593munu_MagDown_full_preselected_iso.root'],
                        'startau':['MC_full_new/Lb_Lc2625taunu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2625taunu_MagDown_full_preselected_iso.root','MC_full_new/Lb_Lc2593taunu_MagUp_full_preselected_iso.root','MC_full_new/Lb_Lc2593taunu_MagDown_full_preselected_iso.root']}
    return startfiles

def GetStartRootFile(startfiles,sample):
    datadir = '/disk/lhcb_data2/RLcMuonic2016/'
    rootfiles = []
    for iF in startfiles[sample]:
        iF = datadir+iF
        rootfiles.append(iF)
    return rootfiles

def WeighSamples(tree,sample):
    if sample=='data':
        return getattr(tree,'sw_sig')
    elif sample=='MISID':
        weight = getattr(tree,'w_recomu_CF')*10*getattr(tree,'sw_sig')
    elif sample=='Combinatorial':
        weight = getattr(tree,'w_MISID')*getattr(tree,'sw_sig')
    else:
        weight = tree.Event_PIDCalibEffWeight
    return weight

def Weight_MultibodyCharm_linear(tree,alpha1):
    weight = 1+2*alpha1*getattr(tree,'w_mbody')
    return weight

def Weight_MultibodyCharm_quadratic(tree,alpha2):
    weight = (1-alpha2)+8*alpha2*getattr(tree,'w_mbody')*getattr(tree,'w_mbody')
    return weight

def GetNormFactor(category,sample,mbody):
    if category=='Kenriched':
        if sample=='2charm' and mbody==False:
            Nc = 365.736
        elif sample=='2charm' and mbody==True:
            Nc = 1574.64
        elif sample == 'Combinatorial':
            Nc = 394.481
        elif sample == 'MISID':
            Nc = 3790.34
        elif sample == 'mu':
            Nc = 13.0116
        elif sample == 'starDs' and mbody==False:
            Nc = 785.547
        elif sample == 'starDs' and mbody==True:
            Nc = 907.587
        elif sample == 'starmu':
            Nc = 408.074
        elif sample=='tau' or sample=='startau':
            Nc=0
    if category=='Isolated':
        if sample=='2charm' and mbody==False:
            Nc = 32666.9
        elif sample=='2charm' and mbody==True:
            Nc = 50892.2
        elif sample == 'Combinatorial':
            Nc = 12777.7
        elif sample == 'MISID':
            Nc =50406
        elif sample == 'mu':
            Nc =801566
        elif sample == 'starDs' and mbody==False:
            Nc = 10
        elif sample == 'starDs' and mbody==True:
            Nc = 5998.76
        elif sample == 'starmu':
            Nc =258808
        elif sample=='tau':
            Nc = 0
        elif sample=='startau':
            Nc=0
    return Nc

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h




def FillHistograms(category):
    hTOT = r.TH1F('hTOT','muPT_TOT',20,0,16E3)
    hTOT_1 = r.TH1F('hTOT_1','muP_TOT',20,0,160E3)
    histos = []
    histos1 = []
    startfiles = GetFileList(category)
    if category=='Kenriched':
        alpha_2charm = [0.481484,0.983569]
        alpha_starDs = [-0.259411,0.535391]
    if category=='Isolated':
        alpha_2charm = [0.0942756,0.41709]
        alpha_starDs = [-2.67839,-0.447195]


    for sample in startfiles.keys():
        print(sample)
        files = GetStartRootFile(startfiles,sample)
        if sample!='2charm' and sample!='starDs' and sample!='data':
            h = r.TH1F('h_'+sample,'muPT_'+sample,20,0,16E3)
            h1 = r.TH1F('h1_'+sample,'muP_'+sample,20,0,160E3)
            h.SetDirectory(0)
            h1.SetDirectory(0)
            for fname in files:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight =WeighSamples(t,sample)
                    if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                        h.Fill(t.mu_PT,weight)
                        h1.Fill(t.mu_P,weight)
            Nc = GetNormFactor(category,sample,False)
            h = ScaleHisto(h,Nc)
            h1 = ScaleHisto(h1,Nc)
            histos.append(h)
            histos1.append(h1)
            hTOT.Add(h)
            hTOT_1.Add(h1)
        elif sample=='2charm' or sample=='starDs':
            h_2b = r.TH1F('h_2b_'+sample,'muPT_'+sample,20,0,16E3)
            h_3b = r.TH1F('h_3b_'+sample,'muPT_'+sample,20,0,16E3)
            h1_2b = r.TH1F('h1_2b_'+sample,'muP_'+sample,20,0,160E3)
            h1_3b = r.TH1F('h1_3b_'+sample,'muP_'+sample,20,0,160E3)
            h_2b.SetDirectory(0)
            h_3b.SetDirectory(0)
            h1_2b.SetDirectory(0)
            h1_3b.SetDirectory(0)
            Nc_2b = GetNormFactor(category,sample,False)
            Nc_3b = GetNormFactor(category,sample,True)
            for fname in files:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight =WeighSamples(t,sample)
                    if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                        if t.twobody==True and t.mbody==False:
                            h_2b.Fill(t.mu_PT,weight)
                            h1_2b.Fill(t.mu_P,weight)
                        if t.twobody==False and t.mbody==True:
                            if sample=='2charm':
                                weight_l = Weight_MultibodyCharm_linear(t,alpha_2charm[0])
                                weight_q = Weight_MultibodyCharm_quadratic(t,alpha_2charm[1])
                            if sample =='starDs':
                                weight_l = Weight_MultibodyCharm_linear(t,alpha_starDs[0])
                                weight_q = Weight_MultibodyCharm_quadratic(t,alpha_starDs[1])
                            h_3b.Fill(t.mu_PT,weight*weight_l*weight_q)
                            h1_3b.Fill(t.mu_P,weight*weight_l*weight_q)
            h_2b = ScaleHisto(h_2b,Nc_2b)
            h1_2b = ScaleHisto(h1_2b,Nc_2b)
            h_3b = ScaleHisto(h_3b,Nc_3b)
            h1_3b = ScaleHisto(h1_3b,Nc_3b)
            histos.append(h_2b)
            histos1.append(h1_2b)
            histos.append(h_3b)
            histos1.append(h1_3b)
            hTOT.Add(h_2b)
            hTOT.Add(h_3b)
            hTOT_1.Add(h1_2b)
            hTOT_1.Add(h1_3b)

        elif sample=='data':
            hPT = r.TH1F('hPT','muPT',20,0,16E3)
            hP = r.TH1F('hP','muP',20,0,160E3)
            hPT.SetDirectory(0)
            hP.SetDirectory(0)
            for fname in files:
                f = r.TFile(fname,'READ')
                t = f.Get('DecayTree')
                for i in range(t.GetEntries()):
                    t.GetEntry(i)
                    weight =WeighSamples(t,sample)
                    if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                        hPT.Fill(t.mu_PT,weight)
                        hP.Fill(t.mu_P,weight)
    return histos, histos1, hTOT, hTOT_1, hPT, hP

if __name__== "__main__":
    #category='Kenriched'
    category='Isolated'
    histos , histos1, hTOT, hTOT1, hPT, hP= FillHistograms(category)

    c = r.TCanvas('c','c',1000,500)
    c.Divide(2,1)
    c.cd(1)
    hPT.SetMarkerStyle(20)
    hPT.SetMarkerSize(1)
    hPT.SetLineColor(r.kBlack)
    hPT.Draw('PE0')
    hTOT.Draw('histo sames')
    c.cd(2)
    hP.SetMarkerStyle(20)
    hP.SetMarkerSize(1)
    hP.SetLineColor(r.kBlack)
    hP.Draw('PE0')
    hTOT1.Draw('histo sames')

    hPTstack = r.THStack('hPTstack','')
    l = r.TLegend(0.75,0.75,0.95,0.9)
    for i,h in enumerate(histos):
        h.SetFillColor(2+i)
        hPTstack.Add(h)
        l.AddEntry(h,h.GetName()[2:],'f')
    c1 = r.TCanvas('c1','c1',500,500)
    hPT.SetMarkerStyle(20)
    hPT.SetMarkerSize(1)
    hPT.SetLineColor(r.kBlack)
    hPT.Draw('PE0')
    hPTstack.Draw('histo sames')
    l.Draw()

    hPstack = r.THStack('hPstack','')
    for i,h in enumerate(histos1):
        h.SetFillColor(2+i)
        hPstack.Add(h)
    c2 = r.TCanvas('c2','c2',500,500)
    hP.SetMarkerStyle(20)
    hP.SetMarkerSize(1)
    hP.SetLineColor(r.kBlack)
    hP.Draw('PE0')
    hPstack.Draw('histo sames')
    l.Draw()
