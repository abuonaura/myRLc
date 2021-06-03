#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT as r

# In[2]:


folder = '/disk/lhcb_data2/RLcMuonic2016/MC_full_trueTrigger/'
mcsamples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds',
             'Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds']
polarities = ['MagUp','MagDown']
variables = {'Lc_PT':['40','0','16000','#Lambda_{c} PT'],'Lc_P':['40','0','200000','#Lambda_{c} P'],
             'FitVar_El_mLc':['40','0','2600','E_{#mu}'],'FitVar_q2_mLc/1.E6':['40','-2','14','q^{2}'],
             'FitVar_Mmiss2_mLc/1.E6':['40','-2','14','M_{miss}^{2}']}


# In[3]:


def PutTogetherPolarityHistos(h):
    h_new = {}
    for mcsample in mcsamples:
        if(h[mcsample]['MagUp']):
            h_new[mcsample] = h[mcsample]['MagUp']
            if(h[mcsample]['MagDown']):
                h_new[mcsample].Add(h[mcsample]['MagDown'])
                h_new[mcsample].SetTitle(mcsample)
                h_new[mcsample].SetDirectory(0)
    return h_new


# In[41]:


def PlotSuperimposed(var,htrue,hemul,Xaxisname, Yaxisname):
    c = r.TCanvas('c_'+var,'',800,1500)
    nsamples = len(mcsamples)
    npads = int(nsamples/2+1)
    c.Divide(2,npads)
    for i,mcsample in enumerate(mcsamples):
        c.cd(i+1)
        htrue[mcsample].SetLineColor(r.kAzure+9)
        htrue[mcsample].GetXaxis().SetTitle(Xaxisname)
        hemul[mcsample].SetLineColor(r.kOrange+2)
        htrue[mcsample].GetYaxis().SetTitle(Yaxisname)
        htrue[mcsample].Draw()
        hemul[mcsample].Draw('sames')
        l = r.TLegend(0.5,0.6,0.8,0.85)
        l.AddEntry(htrue[mcsample],'True Trigger','l')
        l.AddEntry(hemul[mcsample],'Emulated Trigger','l')
        l.Draw('sames')
    return c


# In[5]:


def CompareL0HadronTOS(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wL0TOSEmulation.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_L0HadronDecision_TOS*(TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450)')
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_L0Hadron_TOS_emulated*(TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450)')
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450')
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# In[6]:


def CompareL0GlobalTIS(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wL0TISEmulation.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lb_L0Global_TIS==1&&TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450')
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lb_L0Global_TIS_emulated==1&&(TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450)')
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450')
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# In[7]:


def CompareL0(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wL0TOSEmulation.root','READ')
            ttrig = ftrig.Get('DecayTree')
            ftrig1 = r.TFile(folder+mcsample+'_'+polarity+'_full_wL0TISEmulation.root','READ')
            ttrig1 = ftrig1.Get('DecayTree')
            startt.AddFriend(ttrig)
            startt.AddFriend(ttrig1)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        '(Lc_L0HadronDecision_TOS==1||Lb_L0Global_TIS==1)&&(TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450)')
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        '(Lc_L0Hadron_TOS_emulated==1||Lb_L0Global_TIS_emulated==1)&&(TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450)')
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'TruthMatch==1&&MoreHLT2Sel==1&&nSPDHits<450')
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# In[8]:


def CompareHLT1OneTrack(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    L0cut = '(Lb_L0Global_TIS==1||Lc_L0HadronDecision_TOS == 1)'
    Truthcut = 'MoreHLT2Sel==1&&TruthMatch==1&&nSPDHits<450'
    cut = Truthcut + '&& '+L0cut
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wHLT1OneTrackEmulation_Df.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_Hlt1TrackMVADecision_TOS==1&&'+cut)
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_HLT1TrackMVA_Emu_EffCorrected_TOS==1&&'+cut)
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',cut)
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# In[9]:


def CompareHLT1TwoTracks(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',int(nbins),int(bmin),int(bmax)) 
                          for polarity in polarities} for mcsample in mcsamples}
    L0cut = '(Lb_L0Global_TIS==1||Lc_L0HadronDecision_TOS == 1)'
    Truthcut = 'MoreHLT2Sel==1&&TruthMatch==1&&nSPDHits<450'
    cut = Truthcut + '&& '+L0cut
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wHLT1TwoTracksEmulation.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_Hlt1TwoTrackMVADecision_TOS==1&&'+cut)
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS==1&&'+cut)
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',cut)
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# In[23]:


def CompareHLT1TwoTrackEmuRatio(var,nbins,bmin,bmax):
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemulcorr = {mcsample: {polarity: r.TH1F('hemulcorr_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    L0cut = '(Lb_L0Global_TIS==1||Lc_L0HadronDecision_TOS == 1)'
    Truthcut = 'MoreHLT2Sel==1&&TruthMatch==1&&nSPDHits<450'
    cut = Truthcut + '&& '+L0cut
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wHLT1TwoTracksEmulation.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_HLT1TwoTrackMVA_Emu_TOS==1&&'+cut)
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_emul_corr_'+name+'('+nbins+','+bmin+','+bmax+')',
                        'Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS==1&&'+cut)
            hemulcorr[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_corr_'+name)
            hemulcorr[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
    
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hemulcorr_1 = PutTogetherPolarityHistos(hemulcorr)
    return hemul_1,hemulcorr_1


# In[24]:


def PlotHLT1TwoTrackEmuRatio(var,h,Xaxisname, Yaxisname):
    c = r.TCanvas('c_'+var,'',800,1500)
    nsamples = len(mcsamples)
    npads = int(nsamples/2+1)
    c.Divide(2,npads)
    for i,mcsample in enumerate(mcsamples):
        c.cd(i+1)
        h[mcsample].SetLineColor(r.kAzure+9)
        h[mcsample].GetXaxis().SetTitle(Xaxisname)
        h[mcsample].GetYaxis().SetTitle(Yaxisname)
        h[mcsample].Draw()
    return c


# In[34]:


def CompareHLT1(var,nbins,bmin,bmax):
    htrue = {mcsample: {polarity: r.TH1F('htrue_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    hemul = {mcsample: {polarity: r.TH1F('hemul_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    hnocut = {mcsample: {polarity: r.TH1F('hnocut_'+var+'_'+mcsample+'_'+polarity,'',40,-2,12) 
                          for polarity in polarities} for mcsample in mcsamples}
    L0cut = '(Lb_L0Global_TIS==1||Lc_L0HadronDecision_TOS == 1)'
    Truthcut = 'MoreHLT2Sel==1&&TruthMatch==1&&nSPDHits<450'
    cut = Truthcut + '&& '+L0cut
    for mcsample in mcsamples:
        for polarity in polarities:
            startf = r.TFile(folder+mcsample+'_'+polarity+'_full.root','READ')
            startt = startf.Get('tupleout/DecayTree')
            preself = r.TFile(folder+mcsample+'_'+polarity+'_full_preselectionVars.root','READ')
            preselt = preself.Get('DecayTree')
            startt.AddFriend(preselt)
            ftrig = r.TFile(folder+mcsample+'_'+polarity+'_full_wHLT1OneTrackEmulation_Df.root','READ')
            ttrig = ftrig.Get('DecayTree')
            startt.AddFriend(ttrig)
            ftrig2 = r.TFile(folder+mcsample+'_'+polarity+'_full_wHLT1TwoTracksEmulation.root','READ')
            ttrig2 = ftrig2.Get('DecayTree')
            startt.AddFriend(ttrig2)
            name = mcsample+'_'+polarity
            startt.Draw(var+'>>'+var+'_true_'+name+'('+nbins+','+bmin+','+bmax+')',
                        '(Lc_Hlt1TrackMVADecision_TOS==1||Lc_Hlt1TwoTrackMVADecision_TOS==1)&&'+cut)
            htrue[mcsample][polarity] = r.gPad.GetPrimitive(var+'_true_'+name)
            startt.Draw(var+'>>'+var+'_emul_'+name+'('+nbins+','+bmin+','+bmax+')',
                        '(Lc_HLT1TrackMVA_Emu_EffCorrected_TOS==1||Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS==1)&&'+cut)
            hemul[mcsample][polarity] = r.gPad.GetPrimitive(var+'_emul_'+name)
            startt.Draw(var+'>>'+var+'_nocut_'+name+'('+nbins+','+bmin+','+bmax+')',cut)
            hnocut[mcsample][polarity] = r.gPad.GetPrimitive(var+'_nocut_'+name)
            htrue[mcsample][polarity].SetDirectory(0)
            hemul[mcsample][polarity].SetDirectory(0)
            hnocut[mcsample][polarity].SetDirectory(0)
    
    htrue_1 = PutTogetherPolarityHistos(htrue)
    hemul_1 = PutTogetherPolarityHistos(hemul)
    hnocut_1 = PutTogetherPolarityHistos(hnocut)
    return htrue_1,hemul_1,hnocut_1


# ## L0_HadronTOS emulation comparison

# In[11]:


h_L0H_true = {}
h_L0H_emul = {}
h_tot = {}
heffL0H_true = {}
heffL0H_emul = {}
for key in variables.keys():
    print(key)
    h_L0H_true[key], h_L0H_emul[key], h_tot[key] = CompareL0HadronTOS(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    heffL0H_true[key] = {}
    heffL0H_emul[key] = {}
    for mcsample in mcsamples:
        heffL0H_true[key][mcsample] = h_L0H_true[key][mcsample]
        heffL0H_true[key][mcsample].Divide(heffL0H_true[key][mcsample],h_tot[key][mcsample])
        heffL0H_emul[key][mcsample] = h_L0H_emul[key][mcsample]
        heffL0H_emul[key][mcsample].Divide(heffL0H_emul[key][mcsample],h_tot[key][mcsample])


# In[42]:


cL0H = {}
for key in variables.keys():
    cL0H[key] = PlotSuperimposed(key+'_L0TOS',heffL0H_true[key],heffL0H_emul[key],variables[key][3],
                                 '#epsilon_{L0TOS/ALL}')


# In[44]:


for key in variables.keys():
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cL0H[key].Draw()
    cL0H[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_L0HTOS.png')
    cL0H[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_L0HTOS.root')


# ## Lb_L0Global_TIS

# In[14]:


h_L0G_true = {}
h_L0G_emul = {}
h_tot = {}
heffL0G_true = {}
heffL0G_emul = {}
for key in variables.keys():
    print(key)
    h_L0G_true[key], h_L0G_emul[key], h_tot[key] = CompareL0GlobalTIS(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    heffL0G_true[key] = {}
    heffL0G_emul[key] = {}
    for mcsample in mcsamples:
        heffL0G_true[key][mcsample] = h_L0G_true[key][mcsample]
        heffL0G_true[key][mcsample].Divide(heffL0G_true[key][mcsample],h_tot[key][mcsample])
        heffL0G_emul[key][mcsample] = h_L0G_emul[key][mcsample]
        heffL0G_emul[key][mcsample].Divide(heffL0G_emul[key][mcsample],h_tot[key][mcsample])


# In[15]:


cL0G = {}
for key in variables.keys():
    cL0G[key] = PlotSuperimposed(key+'_L0TIS',heffL0G_true[key],heffL0G_emul[key],variables[key][3],
                                 '#epsilon_{L0TIS/ALL}')


# In[16]:


for key in variables.keys():
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cL0G[key].Draw()
    cL0G[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_L0GTIS.png')


# ## Lb_L0Global_TIS || Lc_L0Hadron_TOS

# In[17]:


h_L0_true = {}
h_L0_emul = {}
h_tot = {}
heffL0_true = {}
heffL0_emul = {}
for key in variables.keys():
    print(key)
    h_L0_true[key], h_L0_emul[key], h_tot[key] = CompareL0(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    heffL0_true[key] = {}
    heffL0_emul[key] = {}
    for mcsample in mcsamples:
        heffL0_true[key][mcsample] = h_L0_true[key][mcsample]
        heffL0_true[key][mcsample].Divide(heffL0_true[key][mcsample],h_tot[key][mcsample])
        heffL0_emul[key][mcsample] = h_L0_emul[key][mcsample]
        heffL0_emul[key][mcsample].Divide(heffL0_emul[key][mcsample],h_tot[key][mcsample])


# In[18]:


cL0 = {}
for key in variables.keys():
    cL0[key] = PlotSuperimposed(key+'_L0',heffL0_true[key],heffL0_emul[key],variables[key][3],
                                 '#epsilon_{L0_{TOS||TIS}/ALL}')
    cL0G[key].Draw()
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cL0G[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_L0.png')


# ## HLT1 One Track

# In[19]:


h_HLT1_1t_true = {}
h_HLT1_1t_emul = {}
h_tot = {}
heffHLT1_1t_true = {}
heffHLT1_1t_emul = {}
for key in variables.keys():
    print(key)
    h_HLT1_1t_true[key], h_HLT1_1t_emul[key], h_tot[key] = CompareHLT1OneTrack(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    heffHLT1_1t_true[key] = {}
    heffHLT1_1t_emul[key] = {}
    for mcsample in mcsamples:
        heffHLT1_1t_true[key][mcsample] = h_HLT1_1t_true[key][mcsample]
        heffHLT1_1t_true[key][mcsample].Divide(heffHLT1_1t_true[key][mcsample],h_tot[key][mcsample])
        heffHLT1_1t_emul[key][mcsample] = h_HLT1_1t_emul[key][mcsample]
        heffHLT1_1t_emul[key][mcsample].Divide(heffHLT1_1t_emul[key][mcsample],h_tot[key][mcsample])


# In[20]:


cHLT11t = {}
for key in variables.keys():
    cHLT11t[key] = PlotSuperimposed(key+'_HLT11t',heffHLT1_1t_true[key],heffHLT1_1t_emul[key],variables[key][3],
                                 '#epsilon_{HLT1_1trk/L0}')
    cHLT11t[key].Draw()
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cHLT11t[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_HLT1OneTrk.png')


# ## HLT1 TwoTrack

# In[21]:


h_HLT1_2t_true = {}
h_HLT1_2t_emul = {}
h_tot = {}
heffHLT1_2t_true = {}
heffHLT1_2t_emul = {}
for key in variables.keys():
    print(key)
    h_HLT1_2t_true[key], h_HLT1_2t_emul[key], h_tot[key] = CompareHLT1TwoTracks(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    heffHLT1_2t_true[key] = {}
    heffHLT1_2t_emul[key] = {}
    for mcsample in mcsamples:
        heffHLT1_2t_true[key][mcsample] = h_HLT1_2t_true[key][mcsample]
        heffHLT1_2t_true[key][mcsample].Divide(heffHLT1_2t_true[key][mcsample],h_tot[key][mcsample])
        heffHLT1_2t_emul[key][mcsample] = h_HLT1_2t_emul[key][mcsample]
        heffHLT1_2t_emul[key][mcsample].Divide(heffHLT1_2t_emul[key][mcsample],h_tot[key][mcsample])


# In[38]:


cHLT12t = {}
for key in variables.keys():
    cHLT12t[key] = PlotSuperimposed(key+'_HLT12t',heffHLT1_2t_true[key],heffHLT1_2t_emul[key],variables[key][3],
                                 '#epsilon_{HLT1_2trks/L0}')
    cHLT12t[key].Draw()
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cHLT12t[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_HLT1TwoTrks.png')


# ## HLT1 Two Tracks Emu Ratio

# In[25]:


h_HLT1_2t_emul = {}
h_HLT1_2t_emulcorr = {}
h_HLT1_2t_ratio = {}
for key in variables.keys():
    print(key)
    h_HLT1_2t_emul[key], h_HLT1_2t_emulcorr[key] = CompareHLT1TwoTrackEmuRatio(key, variables[key][0],variables[key][1],
                                                                      variables[key][2])
    h_HLT1_2t_ratio[key] = {}
    for mcsample in mcsamples:
        h_HLT1_2t_ratio[key][mcsample] = h_HLT1_2t_emul[key][mcsample]
        h_HLT1_2t_ratio[key][mcsample].Divide(h_HLT1_2t_ratio[key][mcsample],h_HLT1_2t_emulcorr[key][mcsample])


# In[28]:


cHLT12tratio = {}
for key in variables.keys():
    cHLT12tratio[key] = PlotHLT1TwoTrackEmuRatio(key+'_ratio',h_HLT1_2t_ratio[key],variables[key][3],
                                 'correction ratio')
    cHLT12tratio[key].Draw()
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cHLT12tratio[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_HLT1TwoTrks_ratio.png')


# ## HLT1 (HLT1_1Track or HLT1_2Tracks)

# In[35]:


h_HLT1_true = {}
h_HLT1_emul = {}
h_tot = {}
heffHLT1_true = {}
heffHLT1_emul = {}
for key in variables.keys():
    print(key)
    h_HLT1_true[key], h_HLT1_emul[key], h_tot[key] = CompareHLT1(key, variables[key][0],variables[key][1],
                                                                 variables[key][2])
    heffHLT1_true[key] = {}
    heffHLT1_emul[key] = {}
    for mcsample in mcsamples:
        heffHLT1_true[key][mcsample] = h_HLT1_true[key][mcsample]
        heffHLT1_true[key][mcsample].Divide(heffHLT1_true[key][mcsample],h_tot[key][mcsample])
        heffHLT1_emul[key][mcsample] = h_HLT1_emul[key][mcsample]
        heffHLT1_emul[key][mcsample].Divide(heffHLT1_emul[key][mcsample],h_tot[key][mcsample])


# In[39]:


cHLT1 = {}
for key in variables.keys():
    cHLT1[key] = PlotSuperimposed(key+'_HLT1',heffHLT1_true[key],heffHLT1_emul[key],variables[key][3],
                                 '#epsilon_{HLT1/L0}')
    cHLT1[key].Draw()
    name = key
    if(key=='FitVar_q2_mLc/1.E6'): name = "FitVar_q2_mLc"
    if(key=='FitVar_Mmiss2_mLc/1.E6'): name = 'FitVar_Mmiss2_mLc'
    cHLT12t[key].SaveAs('plots/TriggerEmulationEfficiencies/'+name+'_HLT1.png')


# In[ ]:




