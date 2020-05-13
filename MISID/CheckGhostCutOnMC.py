import ROOT as r
import sys, os

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC_full/'

if __name__ == "__main__":
    #cuts = [0.3]
    cuts = [0.1,0.15,0.2,0.25,0.3,0.35,0.4]
    polarities = ['MagUp','MagDown']
    #polarities = ['MagUp']
    #types = ['iso','Kenr']
    types = ['iso']
    
    h_El_all = r.TH1F('h_El_all','E*_{#mu}; (MeV);',10, 0., 2600.) 
    h_q2_all = r.TH1F('h_q2_all',"q^{2}", 4, -2.E6, 14.E6)
    h_Mmiss2_all = r.TH1F('h_Mmiss2_all',"M^{2}_{miss}", 10, -2.E6, 14.E6)

    h_El_all_nog = r.TH1F('h_El_all_nog','E*_{#mu}; (MeV);',10, 0., 2600.) 
    h_q2_all_nog = r.TH1F('h_q2_all_nog',"q^{2}", 4, -2.E6, 14.E6)
    h_Mmiss2_all_nog = r.TH1F('h_Mmiss2_all_nog',"M^{2}_{miss}", 10, -2.E6, 14.E6)

    h_muProbNNvsMughost_all = r.TH2F('h_muProbNNvsMughost_all',';muProbNN;muGhost;',100,0,1,25,0,0.5)
    h_muProbNNvsMuP_all = r.TH2F('h_muProbNNvsMuP_all',';muProbNN;muP;',100,0,1,100,0,2.E5)

    h_muP_all = r.TH1F('h_muP_all',';muP(MeV);',25,0,2.E5)
    histos_muP=[]
    histos_muPT=[]
    histos_ratio_P=[]
    for i in range(len(cuts)):
        histos_muP.append(r.TH1F('h_MuP_cut_all_'+str(i),';muP;',25,0,2.E5))
        histos_muPT.append(r.TH1F('h_MuPT_cut_all_'+str(i),';muPT;',30,0,3.E4))
        histos_ratio_P.append(r.TH1F('h_MuP_ratio_all_'+str(i),';muPT;',30,0,3.E4))

    for tt in types:
        print('-------------------------------------------')
        print('         Analysing '+tt+' sample')
        print('-------------------------------------------')
        n=0
        
        
        for cut in cuts:
            n+=1
            if n==1:
                h_muProbNNvsMughost = {}
                h_muProbNNvsMuP = {}
                h_muP={}

            h_MuP_cut = {}
            h_MuPT_cut = {}
            histos_muP[n-1].Reset()
            histos_muPT[n-1].Reset()
            histos_ratio_P[n-1].Reset()

            h_El ={}
            h_q2 ={} 
            h_Mmiss2 ={}
    
            h_El_nog = {}
            h_q2_nog = {}
            h_Mmiss2_nog = {}
        
            h_El_all.Reset()
            h_q2_all.Reset()
            h_Mmiss2_all.Reset()
            
            h_El_all_nog.Reset()
            h_q2_all_nog.Reset()
            h_Mmiss2_all_nog.Reset()
    
            print('-------------------------------------------')
            print('Requiring mu_ProbNNghost < ' + str(cut))
            for polarity in polarities:
                f = r.TFile(filedir + 'Lb_Lcmunu_'+polarity+'_full_preselected_'+tt+'.root')
                t = f.Get('DecayTree')
                df0 = r.RDataFrame(t)
                df0 = df0.Filter('FitVar_El_mLc>0&&FitVar_El_mLc<2600 && FitVar_q2_mLc>-2E6 && FitVar_q2_mLc<14.E6 && FitVar_Mmiss2_mLc>-2.E6 && FitVar_Mmiss2_mLc<14.E6')
                if n==1:
                    h_muP[polarity] = df0.Histo1D(('h_muP_'+polarity,';muP;',25,0,2.E5),"mu_P","Event_PIDCalibEffWeight")
                    h_muProbNNvsMughost[polarity] = df0.Histo2D(("h_muProbNNvsMughost_"+polarity,";mu_ProbNNghost;mu_ghost",100,0,1,25,0,0.5),"mu_ProbNNghost","mu_ghost")
                    h_muProbNNvsMughost_all.Add(h_muProbNNvsMughost[polarity].GetPtr())
                    h_muProbNNvsMuP[polarity] = df0.Histo2D(("h_muProbNNvsMuP_"+polarity,";mu_ProbNNghost;mu_P",100,0,1,100,0,2.E5),"mu_ProbNNghost","mu_P","Event_PIDCalibEffWeight")
                    h_muProbNNvsMuP_all.Add(h_muProbNNvsMuP[polarity].GetPtr())
                    h_muP_all.Add(h_muP[polarity].GetPtr())

                h_El[polarity] = df0.Histo1D(("h_El_"+polarity, "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","Event_PIDCalibEffWeight")
                #h_El_1 = df0.Histo1D(("h_El_1", "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","Event_PIDCalibEffWeight")
                h_q2[polarity] = df0.Histo1D(("h_q2_"+polarity, "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","Event_PIDCalibEffWeight");
                h_Mmiss2[polarity] = df0.Histo1D(("h_Mmiss2_"+polarity, "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","Event_PIDCalibEffWeight");
                h_El_all.Add(h_El[polarity].GetPtr())
                h_q2_all.Add(h_q2[polarity].GetPtr())
                h_Mmiss2_all.Add(h_Mmiss2[polarity].GetPtr())

                df1 = df0.Filter('mu_ProbNNghost<'+str(cut))
                h_El_nog[polarity] = df1.Histo1D(("h_El_nog_"+polarity, "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","Event_PIDCalibEffWeight")
                h_q2_nog[polarity] = df1.Histo1D(("h_q2_nog_"+polarity, "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","Event_PIDCalibEffWeight");
                h_Mmiss2_nog[polarity] = df1.Histo1D(("h_Mmiss2_nog_"+polarity, "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","Event_PIDCalibEffWeight");

                h_El_all_nog.Add(h_El_nog[polarity].GetPtr())
                h_q2_all_nog.Add(h_q2_nog[polarity].GetPtr())
                h_Mmiss2_all_nog.Add(h_Mmiss2_nog[polarity].GetPtr())

                df2 = df0.Filter('mu_ProbNNghost>'+str(cut))
                h_MuP_cut[polarity] = df2.Histo1D(('h_MuP_cut_'+polarity,';muP;',25,0,2.E5),"mu_P","Event_PIDCalibEffWeight")
                histos_muP[n-1].Add(h_MuP_cut[polarity].GetPtr())
        
                h_MuPT_cut[polarity] = df2.Histo1D(('h_MuPT_cut_'+polarity,';muPT;',30,0,3.E4),"mu_PT","Event_PIDCalibEffWeight")
                histos_muPT[n-1].Add(h_MuPT_cut[polarity].GetPtr())
        
            print('Number of events: ', h_El_all.Integral())
            print('Number of events (with ProbNNghost<'+str(cut)+'): ', h_El_all_nog.Integral())
            print('Percentage of events passing the cut {:.3f}'.format(h_El_all_nog.Integral()/ h_El_all.Integral()))
            '''        
            if n==1:
                c = r.TCanvas('c','c',500,500)
                h_muProbNNvsMughost_all.Draw()
                c.SaveAs('plotsMC/muProbNNvsMuGhost_'+tt+'.root')

                c1 = r.TCanvas('c1','c1',500,500)
                h_muProbNNvsMuP_all.Draw()
                c.SaveAs('plotsMC/muProbNNvsMuP_'+tt+'.root')

            c2 = r.TCanvas('c2','c2',1500,500)
            c2.Divide(3,1)
            c2.cd(1)
            h_El_all.Draw('hist')
            h_El_all_nog.SetLineColor(r.kRed)
            h_El_all_nog.Draw('hist same')
            c2.cd(2)
            h_q2_all.Draw("hist")
            h_q2_all_nog.SetLineColor(r.kRed)
            h_q2_all_nog.Draw('hist same')
            c2.cd(3)
            h_Mmiss2_all.Draw("hist")
            h_Mmiss2_all_nog.SetLineColor(r.kRed)
            h_Mmiss2_all_nog.Draw('hist same')
            c2.SaveAs('plotsMC/FitVars_'+tt+'_ProbNNghost_'+str(cut)+'.png')
        '''
        colors = [r.kRed, r.kBlack,r.kGreen+3, r.kBlue-4, r.kViolet, r.kCyan+1, r.kOrange-3]
        c3 = r.TCanvas('c3','c3',500,500)
        leg = r.TLegend(.5,.5,.87,.87)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        c3.SetLogy()
        for i in range(len(histos_muP)):
            histos_muP[i].SetLineColor(colors[i])
            histos_muP[i].Draw('hist same')
            leg.AddEntry(histos_muP[i],'Cut at '+str(cuts[i]),'L')
        leg.Draw()
        c3.SaveAs('plotsMC/MuP_'+tt+'.root')
        
        c4 = r.TCanvas('c4','c4',500,500)
        c4.SetLogy()
        for i in range(len(histos_muPT)):
            histos_muPT[i].SetLineColor(colors[i])
            histos_muPT[i].Draw('hist same')
        leg.Draw()
        c4.SaveAs('plotsMC/MuPT_'+tt+'.root')

        c5 = r.TCanvas('c5','c5',500,500)
        c5.SetLogy()
        h_muP_all.Draw('hist')

        c6 = r.TCanvas('c6','c6',500,500)
        c6.SetLogy()
        for i in range(len(histos_muP)):
            histos_ratio_P[i] = histos_muP[i]
            histos_ratio_P[i].Divide(h_muP_all)
            histos_ratio_P[i].SetLineColor(colors[i])
            histos_ratio_P[i].Draw('hist same')
        leg.Draw()
        c6.SaveAs('plotsMC/RatioP_'+tt+'.root')
