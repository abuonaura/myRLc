import ROOT as r

fdir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_TrueIsoInfo/'
hdir = '/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/'

particles = ['K','Pi']
hcuts = ['IsMuon==0 && DLLK>4.0', 'IsMuon==0 && DLLK>4.0 && DLLp-DLLK<0', 'IsMuon==0 && DLLK>4.0 && DLLp-DLLK>=0']
folder = '/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/'
polarities = ['MagUp','MagDown']
BDTcut = 0.35
BDTcut2 = 0.35
BDTcut3 = 0.35

func_isolation = '''
double GetMuCharge(double muID)
{
    return -muID/13;
}

int GetProbabilityK_0( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";
        
    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_K_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("K_IsMuon==0 && DLLK>4.0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);

    fh->Close();
    //delete fh;
    //delete h;
    return prob;
}

int GetProbabilityK_1( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";
        
    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_K_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("K_IsMuon==0 && DLLK>4.0 && DLLp-DLLK<0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);
    fh->Close();
    return prob;
}

int GetProbabilityK_2( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";
        
    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_K_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("K_IsMuon==0 && DLLK>4.0 && DLLp-DLLK>=0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);
    fh->Close();
    return prob;
}

int GetProbabilityPi_0( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";

    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_Pi_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("Pi_IsMuon==0 && DLLK>4.0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);
    fh->Close();
    return prob;
}

int GetProbabilityPi_1( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";
        
    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_Pi_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("Pi_IsMuon==0 && DLLK>4.0 && DLLp-DLLK<0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);
    fh->Close();
    return prob;
}

int GetProbabilityPi_2( double nTracks, double Lb_ISOLATION_P, double Lb_ISOLATION_eta, int polarity)
{
    string pol;
    if (polarity==1)
        pol="MagUp";
    if (polarity==-1)
        pol = "MagDown";
        
    TFile *fh = new TFile(("/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/PerfHists_Pi_Turbo16_"+pol+"_Brunel_P_Brunel_ETA_nTracks_Brunel.root").c_str());
    TH3D* h = (TH3D*)fh->Get("Pi_IsMuon==0 && DLLK>4.0 && DLLp-DLLK>=0_All");
    int nbin = h->FindFixBin(Lb_ISOLATION_P, Lb_ISOLATION_eta, nTracks);
    double prob=0;
    if(nbin)
        prob = h->GetBinContent(nbin);
    fh->Close();
    return prob;
}

double GetMass(double Lb_ISOLATION_PIDK, double Lb_ISOLATION_CHARGE, double muCharge, double PIDdiff)
{
    double m_pi = 139.57018;
    double m_K = 493.677;
    double m_p = 938.27208;
    double m = m_pi;
    if (Lb_ISOLATION_PIDK>4.)
    {
        if (Lb_ISOLATION_CHARGE==-muCharge || (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff<=0.))
            m=m_K;
        else if (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff>0)
            m=m_p;
    }
    return m;
}

bool CheckIfK(double Lb_ISOLATION_BDT, double Lb_ISOLATION_probK0, double Lb_ISOLATION_probK1, double Lb_ISOLATION_probK2, double NNghost, double Lb_ISOLATION_CHARGE, double Type, double muCharge, double BDTcut)
{
    //TO have always same seed
    TRandom3 *s = new TRandom3();
    double rndm = s->Uniform(0,1);
    bool isK=0;
    if(Lb_ISOLATION_BDT>BDTcut)
    {
        if(Type==3 && NNghost<0.2)
        {
            if(Lb_ISOLATION_CHARGE==-muCharge)
            {
                if(Lb_ISOLATION_probK0>rndm)
                    isK=1;
            }
            else if (Lb_ISOLATION_CHARGE==muCharge)
            {
                if(Lb_ISOLATION_probK1>rndm)
                    isK=1;
            }
        }
    }
    delete s;
    return isK;
}

'''

r.gInterpreter.Declare(func_isolation)

for polarity in polarities:
    fname = fdir+'Lb_LcDs_'+polarity+'_full.root'
    f = r.TFile(fname)
    f_p = r.TFile(fname[0:-5]+'_preselectionVars.root')
    t = f.Get('tupleout/DecayTree')
    t_p = f_p.Get('DecayTree')
    t.AddFriend(t_p)
    df0 = r.RDataFrame(t)

    h = {}
    for particle in particles:
        fhname = folder+'PerfHists_'+particle+'_Turbo16_'+polarity+'_Brunel_P_Brunel_ETA_nTracks_Brunel.root'
        fh = r.TFile(fname,'READ')
        h[particle] = {}
        for hcut in hcuts:
            h[particle][hcut] = fh.Get(particle+'_'+hcut+'_All')

    if polarity =="MagUp":
        df0 = df0.Define('SamplePolarity', str(1))
    if polarity =="MagDown":
        df0 = df0.Define('SamplePolarity',str(-1))
    df0 = df0.Define('Lb_ISOLATION_P',"TMath::Sqrt(Lb_ISOLATION_PX*Lb_ISOLATION_PX+Lb_ISOLATION_PY*Lb_ISOLATION_PY+Lb_ISOLATION_PZ*Lb_ISOLATION_PZ)")
    df0 = df0.Define('Lb_ISOLATION_P2',"TMath::Sqrt(Lb_ISOLATION_PX2*Lb_ISOLATION_PX2+Lb_ISOLATION_PY2*Lb_ISOLATION_PY2+Lb_ISOLATION_PZ2*Lb_ISOLATION_PZ2)")
    df0 = df0.Define('Lb_ISOLATION_P3',"TMath::Sqrt(Lb_ISOLATION_PX3*Lb_ISOLATION_PX3+Lb_ISOLATION_PY3*Lb_ISOLATION_PY3+Lb_ISOLATION_PZ3*Lb_ISOLATION_PZ3)")
    df0 = df0.Define("muCharge","GetMuCharge(mu_ID)")
    df0 = df0.Define('Lb_ISOLATION_eta',"0.5*TMath::Log((Lb_ISOLATION_P+Lb_ISOLATION_PZ)/(Lb_ISOLATION_P-Lb_ISOLATION_PZ))")
    df0 = df0.Define('Lb_ISOLATION_eta2',"0.5*TMath::Log((Lb_ISOLATION_P2+Lb_ISOLATION_PZ2)/(Lb_ISOLATION_P2-Lb_ISOLATION_PZ2))")
    df0 = df0.Define('Lb_ISOLATION_eta3',"0.5*TMath::Log((Lb_ISOLATION_P3+Lb_ISOLATION_PZ3)/(Lb_ISOLATION_P3-Lb_ISOLATION_PZ3))")
    df0 = df0.Define('Lb_ISOLATION_probK0',"GetProbabilityK_0(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK0_2',"GetProbabilityK_0(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK0_3',"GetProbabilityK_0(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK1',"GetProbabilityK_1(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK1_2',"GetProbabilityK_1(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK1_3',"GetProbabilityK_1(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK2',"GetProbabilityK_2(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK2_2',"GetProbabilityK_2(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probK2_3',"GetProbabilityK_2(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi0',"GetProbabilityPi_0(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi0_2',"GetProbabilityPi_0(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi0_3',"GetProbabilityPi_0(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi1',"GetProbabilityPi_1(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi1_2',"GetProbabilityPi_1(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi1_3',"GetProbabilityPi_1(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi2',"GetProbabilityPi_2(nTracks,Lb_ISOLATION_P,Lb_ISOLATION_eta,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi2_2',"GetProbabilityPi_2(nTracks,Lb_ISOLATION_P2,Lb_ISOLATION_eta2,SamplePolarity)")
    df0 = df0.Define('Lb_ISOLATION_probPi2_3',"GetProbabilityPi_2(nTracks,Lb_ISOLATION_P3,Lb_ISOLATION_eta3,SamplePolarity)")

    df0 = df0.Define("BDTcut",str(BDTcut))
    df0 = df0.Define("BDTcut2",str(BDTcut2))
    df0 = df0.Define("BDTcut3",str(BDTcut3))
    df1 = df0.Filter('FinalSel==1&&(Lb_ISOLATION_BDT>'+str(BDTcut)+'||Lb_ISOLATION_BDT2>'+str(BDTcut2)+'||Lb_ISOLATION_BDT3>'+str(BDTcut3)+')')

    #h = df1.Histo1D(('h_ISOLATION_BDT','',50,-1,1),'Lb_ISOLATION_BDT')
    #h = df1.Histo1D(('h_ISOLATION_BDT2','',50,-1,1),'Lb_ISOLATION_BDT2')
    #h = df1.Histo1D(('h_ISOLATION_BDT3','',50,-1,1),'Lb_ISOLATION_BDT3')
    #h.Draw()
    df1 = df1.Define('isK','CheckIfK(Lb_ISOLATION_BDT, Lb_ISOLATION_probK0, Lb_ISOLATION_probK1, Lb_ISOLATION_probK2, Lb_ISOLATION_NNghost, Lb_ISOLATION_CHARGE, Lb_ISOLATION_Type, muCharge, BDTcut)')
    df1 = df1.Define('isK2','CheckIfK(Lb_ISOLATION_BDT2, Lb_ISOLATION_probK0_2, Lb_ISOLATION_probK1_2, Lb_ISOLATION_probK2_2, Lb_ISOLATION_NNghost2, Lb_ISOLATION_CHARGE2, Lb_ISOLATION_Type2, muCharge, BDTcut2)')
    df1 = df1.Define('isK3','CheckIfK(Lb_ISOLATION_BDT3, Lb_ISOLATION_probK0_3, Lb_ISOLATION_probK1_3, Lb_ISOLATION_probK2_3, Lb_ISOLATION_NNghost3, Lb_ISOLATION_CHARGE3, Lb_ISOLATION_Type3, muCharge, BDTcut3)')

    print("Writing the new tree ...")
    new_tree = r.std.vector('string')()
    new_tree.push_back('isK')
    new_tree.push_back('isK2')
    new_tree.push_back('isK3')
    df1.Snapshot("DecayTree","prova.root",new_tree)
    print("Tree written.")
