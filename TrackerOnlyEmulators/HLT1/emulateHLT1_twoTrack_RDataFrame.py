import ROOT as r
import math as m 
from array import array
import sys


func_code = '''
bool GECDecision2(double nVeloClusters, double nITClusters, double nOTClusters, int nSPDHits)
{
    if (nVeloClusters>=6000 || nVeloClusters<=50 || nITClusters>=3000 || nITClusters<=50 || nOTClusters>=15000 || nOTClusters<=50 || nSPDHits >450) 
    {
        return false;
    }
    return true;
}


bool TrackReconstructedDecision2()
{
    //---Get the Track Reconstruction Efficiency correction factor
    double Correction = 1./1.035; //This is a number I get from LHCb-PUB-2015-024
    TRandom3 *s = new TRandom3(0);
    double num = s->Uniform(0.,1.);
    delete s;
    if (num > Correction) {return false;}
    return true;
}

bool TrackInputDecision2(double PT, double P, double TRCHI2, double NDOF, double IP, double GHOSTPROB, double nTTHits)
{
    if (PT<=600 || P<=5000 || TRCHI2/NDOF >= 2.5 || GHOSTPROB>=0.2 || IP<=4. || nTTHits<3.)
    {
        return false;
    }
    else
        return true;
}

bool HLT1TwoTrackMVADecision(double VCHI2, double ETA, double BPVMCORR, double BPVDIRA, double MVA, double SUMPT, double DOCA, double VDCHI2)
{
    //Requirements on track pair before vertexing
    if (SUMPT    <= 2000 || DOCA<=0. || DOCA >= 10.)          {return false;}

    //Requirements on the track pair combination
    if (VCHI2 >= 10. || VCHI2 <= 0 || BPVMCORR <= 1000. || BPVMCORR >= 1000000000.0 || ETA <=2. || ETA>=5. || BPVDIRA  <= 0.)           {return false;}

    // Remove dummmy values from the input to the MVA
    if(VDCHI2 <=0 || SUMPT <=0)  {return false;}

    if (MVA      <= 0.95)         {return false;}
    //std::cout << "Returned true" << std::endl;

    return true;
}

bool emulateHLT1TwoTrackMVA(bool isGECPassed, bool isTrackPassed_K, bool isTrackPassed_p, bool isTrackPassed_pi, 
                            bool isTrackReco_p, bool isTrackReco_K, bool isTrackReco_pi, bool trackRecoCorrection,
                            bool HLT1_2TrackComb_1_2, bool HLT1_2TrackComb_1_3, bool HLT1_2TrackComb_2_3)
{
    //p=1, K=2, pi=3    
    if(isGECPassed==false) return false;
    
    //Check if combination of tracks passes single track cut requirements
    bool isTrackPassed_1_2 = isTrackPassed_p && isTrackPassed_K;
    bool isTrackPassed_1_3 = isTrackPassed_p && isTrackPassed_pi;
    bool isTrackPassed_2_3 = isTrackPassed_K && isTrackPassed_pi;

    //Check if combination of tracks are reconstructed
    bool isTrackReco_1_2, isTrackReco_1_3, isTrackReco_2_3;

    if (trackRecoCorrection==true)
    {
        isTrackReco_1_2 = isTrackReco_p && isTrackReco_K;
        isTrackReco_1_3 = isTrackReco_p && isTrackReco_pi;
        isTrackReco_2_3 = isTrackReco_K && isTrackReco_pi;
    }
    if (trackRecoCorrection==false)
    {
        isTrackReco_1_2 = true;
        isTrackReco_1_3 = true;
        isTrackReco_2_3 = true;
    }

    if(isTrackPassed_1_2 && isTrackReco_1_2 && HLT1_2TrackComb_1_2)
        return true;
    
    if(isTrackPassed_1_3 && isTrackReco_1_3 && HLT1_2TrackComb_1_3)
        return true;

    if(isTrackPassed_2_3 && isTrackReco_2_3 && HLT1_2TrackComb_2_3)
        return true;

    return false;
}

'''


def main(MC_fileName = "Lb_Lcmunu_MagUp.root"):
    F=r.TFile(MC_fileName)
    T=F.Get("tupleout/DecayTree")

    df0 = r.RDataFrame(T)

    r.gInterpreter.Declare(func_code)
    df1 = df0.Define("isGECPassed","GECDecision2(nVeloClusters, nITClusters,nOTClusters, nSPDHits)")
    df1 = df1.Define('isTrackPassed_p',
                     'TrackInputDecision2(Lb_PT_DAU_1, Lb_P_DAU_1, Lb_TRACK_CHI2_DAU_1,Lb_TRACK_NDOF_DAU_1, Lb_IPCHI2_OWNPV_DAU_1, Lb_TRACK_GHOSTPROB_DAU_1, p_TRACK_nTTHits)') #cut in Lb_IPCHI2_OWNPV not in Iaros?
    df1 = df1.Define('isTrackPassed_K',
                     'TrackInputDecision2(Lb_PT_DAU_2, Lb_P_DAU_2, Lb_TRACK_CHI2_DAU_2, Lb_TRACK_NDOF_DAU_2, Lb_IPCHI2_OWNPV_DAU_2, Lb_TRACK_GHOSTPROB_DAU_2, K_TRACK_nTTHits)')
    df1 = df1.Define('isTrackPassed_pi',
                     'TrackInputDecision2(Lb_PT_DAU_3, Lb_P_DAU_3, Lb_TRACK_CHI2_DAU_3,Lb_TRACK_NDOF_DAU_3, Lb_IPCHI2_OWNPV_DAU_3, Lb_TRACK_GHOSTPROB_DAU_3, pi_TRACK_nTTHits)')

    df1 = df1.Define('isTrackReco_p','TrackReconstructedDecision2()')
    df1 = df1.Define('isTrackReco_K','TrackReconstructedDecision2()')
    df1 = df1.Define('isTrackReco_pi','TrackReconstructedDecision2()')

    df1 = df1.Define('SUMPT_1_2',"TMath::Sqrt( (p_PX + K_PX)*(p_PX + K_PX) + (p_PY + K_PY)*(p_PY + K_PY) )")
    df1 = df1.Define('SUMPT_1_3',"TMath::Sqrt( (p_PX + pi_PX)*(p_PX + pi_PX) + (p_PY + pi_PY)*(p_PY + pi_PY) )")
    df1 = df1.Define('SUMPT_2_3',"TMath::Sqrt( (K_PX + pi_PX)*(K_PX + pi_PX) + (K_PY + pi_PY)*(K_PY + pi_PY) )")

    df1 = df1.Define('HLT1_2TrackComb_1_2',"HLT1TwoTrackMVADecision(Lb_VERTEX_CHI2_COMB_1_2, Lb_ETA_COMB_1_2, Lb_MCORR_OWNPV_COMB_1_2, Lb_DIRA_OWNPV_COMB_1_2, Lb_HLt1TwoTrackMVAEmulations_1_2, SUMPT_1_2, Lb_DOCA_COMB_1_2, Lb_VDCHI2_OWNPV_COMB_1_2)")
    df1= df1.Define('HLT1_2TrackComb_1_3','HLT1TwoTrackMVADecision(Lb_VERTEX_CHI2_COMB_1_3, Lb_ETA_COMB_1_3, Lb_MCORR_OWNPV_COMB_1_3, Lb_DIRA_OWNPV_COMB_1_3, Lb_HLt1TwoTrackMVAEmulations_1_3, SUMPT_1_3, Lb_DOCA_COMB_1_3, Lb_VDCHI2_OWNPV_COMB_1_3)')
    df1= df1.Define('HLT1_2TrackComb_2_3','HLT1TwoTrackMVADecision(Lb_VERTEX_CHI2_COMB_2_3, Lb_ETA_COMB_2_3, Lb_MCORR_OWNPV_COMB_2_3, Lb_DIRA_OWNPV_COMB_2_3, Lb_HLt1TwoTrackMVAEmulations_2_3, SUMPT_2_3, Lb_DOCA_COMB_2_3, Lb_VDCHI2_OWNPV_COMB_2_3)')

    '''
    df1 = df1.Define('Lc_HLT1TwoTrackMVA_Emu_TOS', 'emulateHLT1TwoTrackMVA(isGECPassed, \
                            isTrackPassed_K, isTrackPassed_p, isTrackPassed_pi,\
                            isTrackReco_p, isTrackReco_K, isTrackReco_pi, false,\
                            Lb_VERTEX_CHI2_COMB_1_2, Lb_VERTEX_CHI2_COMB_1_3, Lb_VERTEX_CHI2_COMB_2_3,\
                            Lb_ETA_COMB_1_2, Lb_ETA_COMB_1_3,  Lb_ETA_COMB_2_3,\
                            Lb_MCORR_OWNPV_COMB_1_2,Lb_MCORR_OWNPV_COMB_1_3, Lb_MCORR_OWNPV_COMB_2_3,\
                            Lb_DIRA_OWNPV_COMB_1_2, Lb_DIRA_OWNPV_COMB_1_3, Lb_DIRA_OWNPV_COMB_2_3,\
                            Lb_DOCA_COMB_1_2, Lb_DOCA_COMB_1_3, Lb_DOCA_COMB_2_3,\
                            Lb_VDCHI2_OWNPV_COMB_1_2, Lb_VDCHI2_OWNPV_COMB_1_3, Lb_VDCHI2_OWNPV_COMB_2_3,\
                            Lb_IPCHI2_OWNPV_COMB_1_2,Lb_IPCHI2_OWNPV_COMB_1_3, Lb_IPCHI2_OWNPV_COMB_2_3,\
                            Lb_HLt1TwoTrackMVAEmulations_1_2, Lb_HLt1TwoTrackMVAEmulations_1_3, Lb_HLt1TwoTrackMVAEmulations_2_3,\
                            p_PX, p_PY,  p_P, p_PT,\
                            pi_PX, pi_PY,pi_P, pi_PT,\
                            K_PX, K_PY, K_P, K_PT)')
    '''
    df2 = df1.Define('Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS', 'emulateHLT1TwoTrackMVA(isGECPassed, \
                            isTrackPassed_K, isTrackPassed_p, isTrackPassed_pi,\
                            isTrackReco_p, isTrackReco_K, isTrackReco_pi, true,\
                            HLT1_2TrackComb_1_2, HLT1_2TrackComb_1_3, HLT1_2TrackComb_2_3)')

    print("Writing the new tree ...")
    new_tree = r.std.vector('string')()
    #new_tree.push_back("Lc_HLT1TwoTrackMVA_Emu_TOS")
    new_tree.push_back("Lc_HLT1TwoTrackMVA_Emu_EffCorr_TOS")
    df2.Snapshot("DecayTree",MC_fileName[:-5]+"_wHLT1TwoTracksEmulation.root",new_tree)
    print("Tree written.")



