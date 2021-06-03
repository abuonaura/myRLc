import ROOT as r

#Get Performance Histo file
def GetPerfHisto():
    f = r.TFile("/disk/lhcb_data2/RLcMuonic2016/HistoPID/Pi2K/PerfHists_Pi_Turbo16_MagDown_P_ETA_nTracks_Brunel.root",'READ')
    perfhistoname = "Pi_IsMuon==0 && DLLK>4.0_All"
    perfhisto = f.Get(perfhistoname)
    perfhisto.SetDirectory(0)
    return perfhisto

def DrawPerfHisto(h):
    c = r.TCanvas('c','',500,500)
    h.Draw()
    c1 = r.TCanvas('c1','c1',1500,500)
    c1.Divide(3,1)
    c1.cd(1)
    h.ProjectionX().Draw('hist')
    c1.cd(2)
    h.ProjectionY().Draw('hist')
    c1.cd(3)
    h.ProjectionZ().Draw('hist')
    return c, c1

def GetAntiIsolatedPions(sample):
    datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
    f = r.TFile(datadir+sample+'_MagUp.root','READ')
    f1 = r.TFile(datadir+sample+'_MagUp_full_preselectionVars.root','READ')
    f2 = r.TFile(datadir+sample+'_MagUp_full_PIDCalib.root','READ')
    #Get Tree + tree with preselection variables
    t = f.Get('tupleout/DecayTree')
    t1 = f1.Get('DecayTree')
    t2 = f2.Get('tupleout/DecayTree')
    t.AddFriend(t1)
    t.AddFriend(t2)
    #Define histos you want to look at
    h_PiP = r.TH1F('h_PiP','Anti iso true #pi momentum; p(MeV/c);',20,0,5E6)
    h_PiEta = r.TH1F('h_PiEta','Anti iso true #pi eta; #eta;',4,1.5,5)
    h_nTracks = r.TH1F('h_nTracks','Number of tracks in the evt; nTracks;',50,0,1000)
    #Loop over entries
    #for i in range(t.GetEntries()):
    npi=0
    for i in range(100000):
        t.GetEntry(i)
        if t.FinalSel!=1 or t.isIsolated==1:
            continue
        nTracks = t.nTracks
        if abs(t.Lb_ISOLATION_TruePID)==211 and t.Lb_ISOLATION_Type==3 and t.Lb_ISOLATION_NNghost<0.2:
            px = t.Lb_ISOLATION_PX
            py = t.Lb_ISOLATION_PY
            pz = t.Lb_ISOLATION_PZ
            p = r.TMath.Sqrt(px*px+py*py+pz*pz)
            eta = 0.5*r.TMath.Log((p+pz)/(p-pz))
            h_PiP.Fill(p)
            h_PiEta.Fill(eta)
            h_nTracks.Fill(nTracks)
            npi+=1

    print(npi)
    
    c2 = r.TCanvas('c2','c2',1500,500)
    c2.Divide(3,1)
    c2.cd(1)
    h_PiP.Draw()
    c2.cd(2)
    h_PiEta.Draw()
    c2.cd(3)
    h_nTracks.Draw()
    return c2

h = GetPerfHisto()
c, c1 = DrawPerfHisto(h)
c2 = GetAntiIsolatedPions('Lb_LcDs')




