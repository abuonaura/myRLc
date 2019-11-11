'''
Author: Annarita Buonaura
Date: 30 December 2018
Purpose: Draw fit templates of data and misid samples
'''

import ROOT as r

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
polarities=['MagUp','MagDown']
particles=['K','Pi']

datafiles = {'MagUp':datadir+'Data/Lb_DataSS_MagUp_reduced_preselected_sw.root', 'MagDown':datadir+'Data/Lb_DataSS_MagDown_reduced_preselected_sw.root'}

fakemufiles = {'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp_reduced_preselected_sw.root', 'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown_reduced_preselected_sw.root'}

misidfiles = {'MagUp':{'K':datadir+'/MISID/SameSign/K_sample_MagUp_sw.root','Pi':datadir+'/MISID/SameSign/Pi_sample_MagUp_sw.root'},'MagDown':{'K':datadir+'/MISID/SameSign/K_sample_MagDown_sw.root','Pi':datadir+'/MISID/SameSign/Pi_sample_MagDown_sw.root'}}

misidfilesCF = {'MagUp':{'K':datadir+'/MISID/SameSign/K_sample_MagUp_sw_withCF.root','Pi':datadir+'/MISID/SameSign/Pi_sample_MagUp_sw_withCF.root'},'MagDown':{'K':datadir+'/MISID/SameSign/K_sample_MagDown_sw_withCF.root','Pi':datadir+'/MISID/SameSign/Pi_sample_MagDown_sw_withCF.root'}}

def GetLuminosity(fname):
    #print(fname)
    f = r.TFile(fname,'read')
    lt = f.Get('LumiTuple')
    luminosity=0
    for evt in range(lt.GetEntries()+1):
        lt.GetEntry(evt)
        luminosity+=lt.IntegratedLuminosity
    return luminosity

lumiData_Up = GetLuminosity(datafiles['MagUp'])
lumiData_Down = GetLuminosity(datafiles['MagDown'])
lumiData = lumiData_Up+lumiData_Down
lumiFakeMu_Up = GetLuminosity(fakemufiles['MagUp'])
lumiFakeMu_Down = GetLuminosity(fakemufiles['MagDown'])
lumiFakeMu = lumiFakeMu_Up+lumiFakeMu_Down

r_lumiUp  = lumiData_Up/lumiFakeMu_Up
r_lumiDown  = lumiData_Down/lumiFakeMu_Down
r_lumi = {'MagUp':r_lumiUp, 'MagDown':r_lumiDown}
r_lumiTOT = lumiData/lumiFakeMu

h_El_data = r.TH1F('h_El_data',';E_{#mu} (MeV);',20,0,2500)
h_q2_data = r.TH1F('h_q2_data',';q^{2} (MeV^{2});',20,-2E6,18E6)
h_Mmiss_data = r.TH1F('h_Mmiss_data',';Mmiss^{2} (MeV^{2});',20,-2E6,18E6)

def GetTemplatesData():
    for polarity in polarities:
        dataf = r.TFile(datafiles[polarity], 'READ')
        datat = dataf.Get('DecayTree')
        for i in range(datat.GetEntries()):
            datat.GetEntry(i)
            h_El_data.Fill(datat.FitVar_El_mLc,datat.sw_sig)
            h_Mmiss_data.Fill(datat.FitVar_Mmiss2_mLc,datat.sw_sig)
            h_q2_data.Fill(datat.FitVar_q2_mLc,datat.sw_sig)
    return h_El_data, h_q2_data, h_Mmiss_data

h_El_muMISID_K = r.TH1F('h_El_muMISID_K',';E_{#mu} (MeV);',20,0,2500)
h_q2_muMISID_K = r.TH1F('h_q2_muMISID_K',';q^{2} (MeV^{2});',20,-2E6,18E6)
h_Mmiss_muMISID_K = r.TH1F('h_Mmiss_muMISID_K',';Mmiss^{2} (MeV^{2});',20,-2E6,18E6)

h_El_muMISID_Pi = r.TH1F('h_El_muMISID_Pi',';E_{#mu} (MeV);',20,0,2500)
h_q2_muMISID_Pi = r.TH1F('h_q2_muMISID_Pi',';q^{2} (MeV^{2});',20,-2E6,18E6)
h_Mmiss_muMISID_Pi = r.TH1F('h_Mmiss_muMISID_Pi',';Mmiss^{2} (MeV^{2});',20,-2E6,18E6)

h_El_muMISID = {'K':h_El_muMISID_K,'Pi':h_El_muMISID_Pi}
h_q2_muMISID = {'K':h_q2_muMISID_K,'Pi':h_q2_muMISID_Pi}
h_Mmiss_muMISID = {'K':h_Mmiss_muMISID_K,'Pi':h_Mmiss_muMISID_Pi}



def GetTemplatesMISID(particle, crossfeed):
    if crossfeed==0: 
        for polarity in polarities:
            misidf = r.TFile(misidfiles[polarity][particle],'READ')
            misidt = misidf.Get('DecayTree')
            for i in range(misidt.GetEntries()):
                misidt.GetEntry(i)
                h_El_muMISID[particle].Fill(misidt.FitVar_El_mLc,misidt.sw_sig_misid * misidt.w_recomu *10)
                h_q2_muMISID[particle].Fill(misidt.FitVar_q2_mLc,misidt.sw_sig_misid * misidt.w_recomu *10)
                h_Mmiss_muMISID[particle].Fill(misidt.FitVar_Mmiss2_mLc,misidt.sw_sig_misid * misidt.w_recomu *10)
    if crossfeed==1:
        for polarity in polarities:
            misidf = r.TFile(misidfilesCF[polarity][particle],'READ')
            misidt = misidf.Get('DecayTree')
            for i in range(misidt.GetEntries()):
                misidt.GetEntry(i)
                h_El_muMISID[particle].Fill(misidt.FitVar_El_mLc,misidt.sw_sig_misid * misidt.w_recomu_CF *10)
                h_q2_muMISID[particle].Fill(misidt.FitVar_q2_mLc,misidt.sw_sig_misid * misidt.w_recomu_CF *10)
                h_Mmiss_muMISID[particle].Fill(misidt.FitVar_Mmiss2_mLc,misidt.sw_sig_misid * misidt.w_recomu_CF *10)
    h_El_muMISID[particle].SetDirectory(0)
    h_q2_muMISID[particle].SetDirectory(0)
    h_El_muMISID[particle].Draw()
    print(h_El_muMISID[particle].Integral())
    h_Mmiss_muMISID[particle].SetDirectory(0)
    return h_El_muMISID[particle], h_q2_muMISID[particle], h_Mmiss_muMISID[particle]

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h


h_El_data, h_q2_data, h_Mmiss_data = GetTemplatesData()
#h_El_data = ScaleHisto(h_El_data,1)
#h_Mmiss_data = ScaleHisto(h_Mmiss_data,1)
#h_q2_data = ScaleHisto(h_q2_data,1)

'''
for particle in particles:
    h_El_muMISID[particle], h_q2_muMISID[particle], h_Mmiss_muMISID[particle] = GetTemplatesMISID(particle,0)
    h_El_muMISID[particle] = ScaleHisto(h_El_muMISID[particle],r_lumiTOT*h_El_muMISID[particle].Integral())
    h_q2_muMISID[particle] = ScaleHisto(h_q2_muMISID[particle],r_lumiTOT*h_q2_muMISID[particle].Integral())
    h_Mmiss_muMISID[particle] = ScaleHisto(h_Mmiss_muMISID[particle],r_lumiTOT*h_Mmiss_muMISID[particle].Integral())

c = r.TCanvas('c','c',1500,500)
c.Divide(3,1)
c.cd(1)
h_Mmiss_data.Draw('hist')
h_Mmiss_data.SetLineWidth(2)
h_Mmiss_data.SetLineColor(r.kBlack)
h_Mmiss_muMISID_K.SetLineWidth(2)
h_Mmiss_muMISID_K.Draw('hist sames')
h_Mmiss_muMISID_Pi.SetLineWidth(2)
h_Mmiss_muMISID_Pi.SetLineColor(r.kRed)
h_Mmiss_muMISID_Pi.Draw('hist sames')
l = r.TLegend(0.1,0.75,0.4,0.9)
l.AddEntry(h_Mmiss_data,"Data","l")
l.AddEntry(h_Mmiss_muMISID_Pi,"Pi_muMISID","l")
l.AddEntry(h_Mmiss_muMISID_K,"K_muMISID","l")
l.Draw()
c.cd(2)
h_El_data.Draw('hist')
h_El_data.SetLineColor(r.kBlack)
h_El_data.SetLineWidth(2)
h_El_muMISID_K.Draw('hist sames')
h_El_muMISID_K.SetLineWidth(2)
h_El_muMISID_Pi.SetLineColor(r.kRed)
h_El_muMISID_Pi.SetLineWidth(2)
h_El_muMISID_Pi.Draw('hist sames')
c.cd(3)
h_q2_data.Draw('hist')
h_q2_data.SetLineColor(r.kBlack)
h_q2_data.SetLineWidth(2)
h_q2_muMISID_K.Draw('hist sames')
h_q2_muMISID_K.SetLineWidth(2)
h_q2_muMISID_Pi.SetLineColor(r.kRed)
h_q2_muMISID_Pi.SetLineWidth(2)
h_q2_muMISID_Pi.Draw('hist sames')
c.SaveAs('plots/FitTemplates_NoCF.png')

'''


for particle in particles:
    h_El_muMISID[particle], h_q2_muMISID[particle], h_Mmiss_muMISID[particle] = GetTemplatesMISID(particle,0)
    h_El_muMISID[particle] = ScaleHisto(h_El_muMISID[particle],r_lumiTOT*h_El_muMISID[particle].Integral())
    h_q2_muMISID[particle] = ScaleHisto(h_q2_muMISID[particle],r_lumiTOT*h_q2_muMISID[particle].Integral())
    h_Mmiss_muMISID[particle] = ScaleHisto(h_Mmiss_muMISID[particle],r_lumiTOT*h_Mmiss_muMISID[particle].Integral())

c1 = r.TCanvas('c','c',1500,500)
c1.Divide(3,1)
c1.cd(1)
h_Mmiss_data.Draw('hist')
h_Mmiss_data.SetLineWidth(2)
h_Mmiss_data.SetLineColor(r.kBlack)
h_Mmiss_muMISID_K.SetLineWidth(2)
h_Mmiss_muMISID_K.Draw('hist sames')
h_Mmiss_muMISID_Pi.SetLineWidth(2)
h_Mmiss_muMISID_Pi.SetLineColor(r.kRed)
h_Mmiss_muMISID_Pi.Draw('hist sames')
l = r.TLegend(0.1,0.75,0.4,0.9)
l.AddEntry(h_Mmiss_data,"Data","l")
l.AddEntry(h_Mmiss_muMISID_Pi,"Pi_muMISID","l")
l.AddEntry(h_Mmiss_muMISID_K,"K_muMISID","l")
l.Draw()
c1.cd(2)
h_El_data.Draw('hist')
h_El_data.SetLineColor(r.kBlack)
h_El_data.SetLineWidth(2)
h_El_muMISID_K.Draw('hist sames')
h_El_muMISID_K.SetLineWidth(2)
h_El_muMISID_Pi.SetLineColor(r.kRed)
h_El_muMISID_Pi.SetLineWidth(2)
h_El_muMISID_Pi.Draw('hist sames')
c1.cd(3)
h_q2_data.Draw('hist')
h_q2_data.SetLineColor(r.kBlack)
h_q2_data.SetLineWidth(2)
h_q2_muMISID_K.Draw('hist sames')
h_q2_muMISID_K.SetLineWidth(2)
h_q2_muMISID_Pi.SetLineColor(r.kRed)
h_q2_muMISID_Pi.SetLineWidth(2)
h_q2_muMISID_Pi.Draw('hist sames')
c1.SaveAs('plots/FitTemplates_WithCF.png')

