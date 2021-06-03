import ROOT as r

polarities = ['MagUp','MagDown']
categories = ['Isolated','Kenriched','Lcpipi']
suffix = {'Isolated':'iso','Kenriched':'Kenr','Lcpipi':'Lcpipi'}
ddir = '/disk/lhcb_data2/RLcMuonic2016/Data/'

h = r.TH3F('h',"qem",4,-2,14,10,0,2600,10,-2,14)
for category in categories:
    h.Reset()
    for polarity in polarities:
        fname = ddir + 'Lb_Data_'+polarity+'_preselected_'+suffix[category]+'_sw.root'
        f = r.TFile(fname,'READ')
        t = f.Get('DecayTree')
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            weight = t.sw_sig
            if (t.FitVar_q2_mLc/1.E6>-2 and t.FitVar_q2_mLc/1.E6<14) and (t.FitVar_Mmiss2_mLc/1.E6>-2 and t.FitVar_Mmiss2_mLc/1.E6<14) and (t.FitVar_El_mLc>0 and t.FitVar_El_mLc<2600):
                h.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,weight)
    print('Category: ', category)
    print('    Signal Yield: ', h.Integral())
