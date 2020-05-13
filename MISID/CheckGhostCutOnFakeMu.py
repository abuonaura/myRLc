import ROOT as r
import sys, os

filedir = '/disk/lhcb_data2/RLcMuonic2016/Data/'

if __name__ == "__main__":
    #cuts = [0.3]
    cuts = [0.1,0.15,0.2,0.25,0.3,0.35,0.4]
    polarities = ['MagUp','MagDown']
    types = ['iso','Kenr']
    
    h_El_all = r.TH1F('h_El_all','E*_{#mu}; (MeV);',10, 0., 2600.) 
    h_q2_all = r.TH1F('h_q2_all',"q^{2}", 4, -2.E6, 14.E6)
    h_Mmiss2_all = r.TH1F('h_Mmiss2_all',"M^{2}_{miss}", 10, -2.E6, 14.E6)

    h_El_all_nog = r.TH1F('h_El_all_nog','E*_{#mu}; (MeV);',10, 0., 2600.) 
    h_q2_all_nog = r.TH1F('h_q2_all_nog',"q^{2}", 4, -2.E6, 14.E6)
    h_Mmiss2_all_nog = r.TH1F('h_Mmiss2_all_nog',"M^{2}_{miss}", 10, -2.E6, 14.E6)

    h_muProbNNvsMughost_all = r.TH2F('h_muProbNNvsMughost_all',';muProbNN;muGhost;',100,0,1,50,0,0.5)
    for tt in types:
        print('Analysing '+tt+' sample')
        n=0
        for cut in cuts:
            n+=1
            if n==1:
                h_muProbNNvsMughost = {}

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
                f = r.TFile(filedir + 'Lb_FakeMu_'+polarity+'_preselected_'+tt+'_sw.root')
                t = f.Get('DecayTree')
                df0 = r.RDataFrame(t)
                df0 = df0.Filter('FitVar_El_mLc>0&&FitVar_El_mLc<2600 && FitVar_q2_mLc>-2E6 && FitVar_q2_mLc<14.E6 && FitVar_Mmiss2_mLc>-2.E6 && FitVar_Mmiss2_mLc<14.E6')
                if n==1:
                    h_muProbNNvsMughost[polarity] = df0.Histo2D(("h_muProbNNvsMughost_"+polarity,";mu_ProbNNghost;mu_ghost",100,0,1,50,0,0.5),"mu_ProbNNghost","mu_ghost")
                    h_muProbNNvsMughost_all.Add(h_muProbNNvsMughost[polarity].GetPtr())

                h_El[polarity] = df0.Histo1D(("h_El_"+polarity, "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","sw_sig")
                #h_El_1 = df0.Histo1D(("h_El_1", "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","sw_sig")
                h_q2[polarity] = df0.Histo1D(("h_q2_"+polarity, "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","sw_sig");
                h_Mmiss2[polarity] = df0.Histo1D(("h_Mmiss2_"+polarity, "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","sw_sig");
                h_El_all.Add(h_El[polarity].GetPtr())
                h_q2_all.Add(h_q2[polarity].GetPtr())
                h_Mmiss2_all.Add(h_Mmiss2[polarity].GetPtr())

                df1 = df0.Filter('mu_ProbNNghost<'+str(cut))
                h_El_nog[polarity] = df1.Histo1D(("h_El_nog_"+polarity, "E*_{#mu}", 10, 0., 2600.),"FitVar_El_mLc","sw_sig")
                h_q2_nog[polarity] = df1.Histo1D(("h_q2_nog_"+polarity, "q^{2}", 4, -2.E6, 14.E6),"FitVar_q2_mLc","sw_sig");
                h_Mmiss2_nog[polarity] = df1.Histo1D(("h_Mmiss2_nog_"+polarity, "M^{2}_{miss}", 10, -2.E6, 14.E6),"FitVar_Mmiss2_mLc","sw_sig");
                
                h_El_all_nog.Add(h_El_nog[polarity].GetPtr())
                h_q2_all_nog.Add(h_q2_nog[polarity].GetPtr())
                h_Mmiss2_all_nog.Add(h_Mmiss2_nog[polarity].GetPtr())
        
        
        
            print('Number of events: ', h_El_all.Integral())
            print('Number of events (with ProbNNghost<'+str(cut)+'): ', h_El_all_nog.Integral())
            print('Percentage of events passing the cut {:.3f}'.format(h_El_all_nog.Integral()/ h_El_all.Integral()))
                
            if n==1:
                c = r.TCanvas('c','c',500,500)
                h_muProbNNvsMughost_all.Draw()
                c.SaveAs('plots/FakeMu_muProbNNvsMuGhost_'+tt+'.png')
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
            c2.SaveAs('plots/FitVars_FakeMu_'+tt+'_ProbNNghost_'+str(cut)+'.png')
