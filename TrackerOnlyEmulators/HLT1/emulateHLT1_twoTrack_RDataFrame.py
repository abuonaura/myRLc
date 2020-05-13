import ROOT as r
import math as m 
from array import array
import sys

func_code =  '''
bool PTDes(double p_PT, double pi_PT,double K_PT)
{
    std::array<double,3> PT = {p_PT,pi_PT,K_PT};
    std::sort(PT.begin(), PT.end());
    if(PT[2]>600 && PT[1]>600) return true;
    return false;
}

bool PDes(double p_P, double pi_P,double K_P)
{
    std::array<double,3> PT = {p_P,pi_P,K_P};
    std::sort(PT.begin(), PT.end());
    if(PT[2]>5000 && PT[1]>5000) return true;
    return false;
}

bool SUMPTDes(int nSPDHits,double nVeloClusters,  double nITClusters, 
                double nOTClusters,
                double p_PX, double p_PY, double p_P, double p_PT, int p_TRACK_nTTHits,
                double pi_PX, double pi_PY, double pi_P, double pi_PT, int pi_TRACK_nTTHits,
                double K_PX, double K_PY, double K_P, double K_PT, int K_TRACK_nTTHits,
                double Lb_IPCHI2_OWNPV_COMB_1_2,double Lb_IPCHI2_OWNPV_COMB_1_3, double Lb_IPCHI2_OWNPV_COMB_2_3,
                double Lb_TRACK_CHI2_DAU_1, double Lb_TRACK_NDOF_DAU_1,
                double Lb_TRACK_CHI2_DAU_2, double Lb_TRACK_NDOF_DAU_2,
                double Lb_TRACK_CHI2_DAU_3, double Lb_TRACK_NDOF_DAU_3,
                double Lb_HLt1TwoTrackMVAEmulations_1_2,double Lb_HLt1TwoTrackMVAEmulations_1_3, double Lb_HLt1TwoTrackMVAEmulations_2_3,
                double Lb_TRACK_GHOSTPROB_DAU_1, double Lb_TRACK_GHOSTPROB_DAU_2, double Lb_TRACK_GHOSTPROB_DAU_3                
                )
{   
    if(nVeloClusters >= 6000 || nITClusters >=3000 || nOTClusters>=15000 || nVeloClusters <= 50 || nITClusters <=50 || nOTClusters<=50 || nSPDHits >450) return false;
    double SUMPT_1_2 = TMath::Sqrt( (p_PX + K_PX)*(p_PX + K_PX) + (p_PY + K_PY)*(p_PY + K_PY) );
    double SUMPT_1_3 = TMath::Sqrt( (p_PX + pi_PX)*(p_PX + pi_PX) + (p_PY + pi_PY)*(p_PY + pi_PY) );
    double SUMPT_2_3 = TMath::Sqrt( (K_PX + pi_PX)*(K_PX + pi_PX) + (K_PY + pi_PY)*(K_PY + pi_PY) );
    
    

    int decision = 0 ;
    if (SUMPT_1_2>2000)
    {
        if (p_P>5000 && K_P>5000)
        {
            if (p_PT>600 && K_PT>600)
            {
                if (Lb_IPCHI2_OWNPV_COMB_1_2 > 4)
                {
                    if ((Lb_TRACK_CHI2_DAU_1/Lb_TRACK_NDOF_DAU_1) < 2.5 && (Lb_TRACK_CHI2_DAU_2/Lb_TRACK_NDOF_DAU_2) < 2.5)
                    {
                        if (Lb_HLt1TwoTrackMVAEmulations_1_2>0.95)
                        {
                            if(Lb_TRACK_GHOSTPROB_DAU_1<0.2 && Lb_TRACK_GHOSTPROB_DAU_2<0.2)
                            {
                                if (p_TRACK_nTTHits>2 && K_TRACK_nTTHits>2)
                                {
                                    decision = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

if (decision ==0 )
{

    if (SUMPT_1_3>2000)
    {
        if (p_P>5000 && pi_P>5000)
        {
            if (p_PT>600 && pi_PT>600)
            {
                if (Lb_IPCHI2_OWNPV_COMB_1_3 > 4)
                {
                    if ((Lb_TRACK_CHI2_DAU_1/Lb_TRACK_NDOF_DAU_1) < 2.5 && (Lb_TRACK_CHI2_DAU_3/Lb_TRACK_NDOF_DAU_3) < 2.5)
                    {
                        if (Lb_HLt1TwoTrackMVAEmulations_1_3>0.95)
                        {
                            if(Lb_TRACK_GHOSTPROB_DAU_1<0.2 && Lb_TRACK_GHOSTPROB_DAU_3<0.2)
                            {
                                if (p_TRACK_nTTHits>2 && pi_TRACK_nTTHits>2)
                                {
                                    decision = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}


if (decision ==0 )
{

    if (SUMPT_2_3>2000)
    {
        if (K_P>5000 && pi_P>5000)
        {
            if (K_PT>600 && pi_PT>600)
            {
                if (Lb_IPCHI2_OWNPV_COMB_2_3 > 4)
                {
                    if ((Lb_TRACK_CHI2_DAU_2/Lb_TRACK_NDOF_DAU_2) < 2.5 && (Lb_TRACK_CHI2_DAU_3/Lb_TRACK_NDOF_DAU_3) < 2.5)
                    {
                        if (Lb_HLt1TwoTrackMVAEmulations_2_3>0.95)
                        {
                            if(Lb_TRACK_GHOSTPROB_DAU_2<0.2 && Lb_TRACK_GHOSTPROB_DAU_3<0.2)
                            {
                                if (K_TRACK_nTTHits>2 && pi_TRACK_nTTHits>2)
                                {
                                    decision = 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}




    if (decision==1) return true;
    else return false;
}

bool TDes(bool P, bool PT, bool SUMPT)
{
    if (P==true && PT==true && SUMPT ==true) return true;
    else return false;
}
'''


def main(MC_fileName = "Lb_Lcmunu_MagUp.root"):
    F=r.TFile(MC_fileName)
    T=F.Get("tupleout/DecayTree")


    r.gInterpreter.Declare(func_code)
    df0 = r.RDataFrame(T)
    df1 = df0.Define("PTDes" ,  "PTDes(p_PT, pi_PT,K_PT)")
    df2 = df1.Define("PDes","PDes(p_P, pi_P, K_P)")
    df3 = df2.Define("SUMPTDes", "SUMPTDes(nSPDHits,nVeloClusters, nITClusters, nOTClusters, \
                p_PX, p_PY, p_P, p_PT, p_TRACK_nTTHits, \
                pi_PX, pi_PY, pi_P, pi_PT, pi_TRACK_nTTHits, \
                K_PX, K_PY, K_P, K_PT, K_TRACK_nTTHits, \
                Lb_IPCHI2_OWNPV_COMB_1_2, Lb_IPCHI2_OWNPV_COMB_1_3, Lb_IPCHI2_OWNPV_COMB_2_3, \
                Lb_TRACK_CHI2_DAU_1, Lb_TRACK_NDOF_DAU_1, \
                Lb_TRACK_CHI2_DAU_2, Lb_TRACK_NDOF_DAU_2, \
                Lb_TRACK_CHI2_DAU_3, Lb_TRACK_NDOF_DAU_3, \
                Lb_HLt1TwoTrackMVAEmulations_1_2, Lb_HLt1TwoTrackMVAEmulations_1_3, Lb_HLt1TwoTrackMVAEmulations_2_3, \
                Lb_TRACK_GHOSTPROB_DAU_1, Lb_TRACK_GHOSTPROB_DAU_2, Lb_TRACK_GHOSTPROB_DAU_3 \
                )")
    df4 = df3.Define("Lc_HLT1TwoTrackMVA_Emu_TOS","TDes(PDes,PTDes,SUMPTDes)")
    print("Writing the new tree ...")
    new_tree = r.std.vector('string')()
    new_tree.push_back("Lc_HLT1TwoTrackMVA_Emu_TOS")
    df4.Snapshot("DecayTree",MC_fileName[:-5]+"_wHLT1TwoTracksEmulation.root",new_tree)
    print("Tree written.")
