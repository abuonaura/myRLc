import ROOT as r
from ROOT import TMath as mt

#Starting variables
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)
ISOBDTcut=0.35
ISOBDT2cut=0.2

datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
samples = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
#samples = ['Lcmunu']
samples = ['LcDs']

func_isolation = '''
double GetMuCharge(double muID)
{
    double muCharge = -1.*muID/13;
    return muCharge;
}

double GetMass(double Lb_ISOLATION_PIDK, double Lb_ISOLATION_CHARGE, double muCharge, double PIDdiff)
{
    double m_pi = 139.57018;
    double m_K = 493.677;
    double m_p = 938.27208;
    double m = m_pi;
    if (Lb_ISOLATION_PIDK>4.)
    {
        if (Lb_ISOLATION_CHARGE==-muCharge || (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff<=0))
            m=m_K;
        else if (Lb_ISOLATION_CHARGE==muCharge&& PIDdiff>0)
            m=m_p;
    }
    return m;
}

bool CheckIfKenriched(double Lb_ISOLATION_BDT, double ISOBDTcut, double Lb_ISOLATION_BDT2, double ISOBDT2cut, double mLc12, double mTOT, double m1, double m2, double Type, double Type2, double NNghost, double NNghost2)
{
    double m_K = 493.677;
    if (Type==3 && NNghost<0.2 && Type2==3 && NNghost2<0.2)
    {
        if(Lb_ISOLATION_BDT>ISOBDTcut && Lb_ISOLATION_BDT2>ISOBDT2cut)
        {
            //cout<<mLc12<<" "<< mTOT<<" "<< m1<<" "<< m2<<endl;
            if (mLc12>2700 && mTOT<5620 && (m1==m_K || m2==m_K))
            {   
                //cout<<"ACCEPTED"<<endl;
                return 1;
            }
        }
    }
    return 0;
}

bool isLcpipi(bool isKenriched, double Lb_ISOLATION_BDT, double ISOBDTcut, double Lb_ISOLATION_BDT2, double ISOBDT2cut, double mLc12, double mTOT, double Type, double Type2, double NNghost, double NNghost2, double charge, double charge2)
{
    if (Type==3 && NNghost<0.2 && Type2==3 && NNghost2<0.2)
    {
        if(Lb_ISOLATION_BDT>ISOBDTcut && Lb_ISOLATION_BDT2>ISOBDT2cut && isKenriched==false && charge==-charge2)
        {
            if(mLc12<2700)
            {
                return true;
                }
            else
                return false;
        }
        else
            return false;
    }
    else
        return false;
}
'''

#interprect c++ code
r.gInterpreter.Declare(func_isolation)

for sample in samples:
    print('--------------- '+sample+' ---------------')
    f = r.TFile(datadir+'Lb_'+sample+'_MagUp.root','READ')
    f1 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_preselectionVars.root')
    f2 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_PIDCalib.root')

    #Get Tree + tree with preselection variables
    t = f.Get('tupleout/DecayTree')
    t1 = f1.Get('DecayTree')
    t2 = f2.Get('tupleout/DecayTree')

    t.AddFriend(t1)
    t.AddFriend(t2)

    #Get dataframe
    df = r.RDataFrame(t)

    #print column names
    '''
    colnames = df.GetColumnNames()
    for col in colnames:
        print(col)
    '''



    df0 = df.Filter("FinalSel==true")

    df1 = df0.Define("E1pi","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX,2)+TMath::Power(Lb_ISOLATION_PY,2)+TMath::Power(Lb_ISOLATION_PZ,2)+"+str(m_pi)+"*"+str(m_pi)+")")
    df1 = df1.Define("E2pi","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX2,2)+TMath::Power(Lb_ISOLATION_PY2,2)+TMath::Power(Lb_ISOLATION_PZ2,2)+"+str(m_pi)+"*"+str(m_pi)+")")
    df1 = df1.Define("ELc","TMath::Sqrt(TMath::Power(Lc_PX,2)+TMath::Power(Lc_PY,2)+TMath::Power(Lc_PZ,2)+"+str(m_Lc)+"*"+str(m_Lc)+")")
    df1 = df1.Define("pLc12_x","Lc_PX+Lb_ISOLATION_PX+Lb_ISOLATION_PX2")
    df1 = df1.Define("pLc12_y","Lc_PY+Lb_ISOLATION_PY+Lb_ISOLATION_PY2")
    df1 = df1.Define("pLc12_z","Lc_PZ+Lb_ISOLATION_PZ+Lb_ISOLATION_PZ2")
    df1 = df1.Define("mLc12","TMath::Sqrt(TMath::Power(ELc+E1pi+E2pi,2)-(TMath::Power(pLc12_x,2)+TMath::Power(pLc12_y,2)+TMath::Power(pLc12_z,2)))")
    df1 = df1.Define("muCharge","GetMuCharge(mu_ID)")
    df1 = df1.Define("PIDdiff","Lb_ISOLATION_PIDp - Lb_ISOLATION_PIDK")
    df1 = df1.Define("PIDdiff2","Lb_ISOLATION_PIDp2 - Lb_ISOLATION_PIDK2")
    df1 = df1.Define("m1","GetMass(Lb_ISOLATION_PIDK, Lb_ISOLATION_CHARGE, muCharge,PIDdiff)")
    df1 = df1.Define("m2","GetMass(Lb_ISOLATION_PIDK2, Lb_ISOLATION_CHARGE2, muCharge,PIDdiff2)")
    df1 = df1.Define("E1","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX,2)+TMath::Power(Lb_ISOLATION_PY,2)+TMath::Power(Lb_ISOLATION_PZ,2)+TMath::Power(m1,2))")
    df1 = df1.Define("E2","TMath::Sqrt(TMath::Power(Lb_ISOLATION_PX2,2)+TMath::Power(Lb_ISOLATION_PY2,2)+TMath::Power(Lb_ISOLATION_PZ2,2)+TMath::Power(m2,2))")
    df1 = df1.Define("Emu","TMath::Sqrt(TMath::Power(mu_PX,2)+TMath::Power(mu_PY,2)+TMath::Power(mu_PZ,2)+"+str(m_mu)+"*"+str(m_mu)+")")
    df1 = df1.Define("Etot","ELc+E1+E2+Emu")
    df1 = df1.Define("pTOT_x","pLc12_x+mu_PX")
    df1 = df1.Define("pTOT_y","pLc12_y+mu_PY")
    df1 = df1.Define("pTOT_z","pLc12_z+mu_PZ")
    df1 = df1.Define("mTOT","TMath::Sqrt(TMath::Power(Etot,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)")
    df2 = df1.Define("SelNew","CheckIfKenriched(Lb_ISOLATION_BDT,"+str(ISOBDTcut)+", Lb_ISOLATION_BDT2,"+str(ISOBDT2cut)+", mLc12, mTOT, m1,  m2, Lb_ISOLATION_Type, Lb_ISOLATION_Type2, Lb_ISOLATION_NNghost, Lb_ISOLATION_NNghost2)")
    #h = df2.Histo1D(('h','',2,0,1),'SelNew')
    #h.Draw()
    '''
    colnames = df1.GetColumnNames()
    for col in colnames:
        print(col)
    '''

    h_mLc12 = df2.Histo1D(('h_mLc12','',100,2000,5000),'mLc12')
    h_E1pi = df2.Histo1D(('h_E1pi','',100,0,10000),'E1pi')
    h_E2pi = df2.Histo1D(('h_E2pi','',100,0,10000),'E2pi')
    h_ELc = df2.Histo1D(('h_ELc','',100,0,1E6),'ELc')
    h_m1 = df2.Histo1D(('h_m1','',100,0,1000),'m1')
    h_m2 = df2.Histo1D(('h_m2','',100,0,1000),'m2')
    h_muCh = df2.Histo1D(('h_muCh','',3,-1,1),'muCharge')
    h_PIDdiff = df2.Histo1D(('h_PIDdiff','',100,-1000,1000),'PIDdiff')
    h_PIDdiff2 = df2.Histo1D(('h_PIDdiff2','',100,-1000,1000),'PIDdiff2')

    c4 = r.TCanvas('c4','c4',500,500)
    h_mLc12.Draw()
    
    c5 = r.TCanvas('c5','c5',1500,500)
    c5.Divide(3,1)
    c5.cd(1)
    h_E1pi.Draw('hist')
    c5.cd(2)
    h_E2pi.Draw('hist')
    c5.cd(3)
    h_ELc.Draw('hist') 

    c6 = r.TCanvas('c6','c6',1000,500)
    c6.Divide(2,1)
    c6.cd(1)
    h_m1.Draw('hist')
    c6.cd(2)
    h_m2.Draw('hist')

    c7 = r.TCanvas('c7','c7',500,500)
    h_muCh.Draw()

    c8 = r.TCanvas('c8','c8',1000,500)
    c8.Divide(2,1)
    c8.cd(1)
    h_PIDdiff.Draw()
    c8.cd(2)
    h_PIDdiff2.Draw()

    entries = df2.Filter("SelNew==true").Count()
    print(entries.GetValue())



