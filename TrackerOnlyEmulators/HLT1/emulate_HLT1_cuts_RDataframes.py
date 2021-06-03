
from array import array
import sys


func_code = '''
bool GECDecision(double nVeloClusters, double nITClusters, double nOTClusters)
{
    if (nVeloClusters>=6000 || nVeloClusters<=50 || nITClusters>=3000 || nITClusters<=50 || nOTClusters>=15000 ||
        nOTClusters<=50) 
    {
        return false;
    }
    return true;
}

bool TrackInputDecision(double PT, double P, double TRCHI2, double GHOSTPROB, double NDOF)
{
    if (PT <= 500 || P  <= 5000 || TRCHI2/NDOF >= 4.0 || GHOSTPROB>999.)     
    {
        return false;
    }
    return true;
}

bool TrackAdditionalDecision(double nTTHits)
{
    if(nTTHits<3.) 
        return false;
    return true;
}

bool TrackMVADecision(double PT, double TRCHI2, double GHOSTPROB, double IPCHI2, double NDOF)
{
    double CHI2 = TRCHI2/NDOF;
    if(PT<500 || CHI2>2.5) {return false;}
    if(GHOSTPROB>=0.2) {return false;}
    if(PT>25000)
    {
        if(IPCHI2>7.4)
            return true;
        else
            return false;
    }
    if(IPCHI2>0 && PT!=500)
    {
        if (TMath::Log(IPCHI2)>1.0/ TMath::Power(PT/1000. - 1.0,2) + 1.1/25000.*(25000. - PT) + TMath::Log(7.4))
            return true;
    }
    return false;
}

bool Inefficiency_Decision(double num, double Correction)
{
    if (num > Correction){return false;}
    return true;
}

bool TrackReconstructedDecision()
{
    //---Get the Track Reconstruction Efficiency correction factor
    double Correction = 1./1.035; //This is a number I get from LHCb-PUB-2015-024
    TRandom3 *s = new TRandom3(0);
    double num = s->Uniform(0.,1.);
    bool isReco =  Inefficiency_Decision(num,Correction);
    delete s;
    if(isReco==0) {return false;}
    return true;
}

bool HLT1TrackMVA_TOS( bool isGECPassed, bool isTrackPassed, bool isAdditional,bool isPassed)
{
    //cout<<"isGECPassed: "<<isGECPassed<<" isTrackPassed: "<<isTrackPassed<<"  isAdditional: "<<isAdditional<<" isPassed:"<< isPassed<<endl;
    if (isPassed && isGECPassed && isTrackPassed && isAdditional)
        return true;
    else
        return false;
}

bool HLT1TrackMVA_TOS_EffCorrected(bool isGECPassed, bool isTrackPassed, bool isAdditional,bool isPassed, bool isTrackReco)
{
    if (isTrackReco && isPassed && isGECPassed && isTrackPassed && isAdditional)
        return true;
    else
        return false;
}

bool EmulationHLT1(bool HLT1TrackMVA_TOS_K,bool HLT1TrackMVA_TOS_p,bool HLT1TrackMVA_TOS_pi)
{
    if(HLT1TrackMVA_TOS_K||HLT1TrackMVA_TOS_p||HLT1TrackMVA_TOS_pi)
        return true;
    else
        return false;
}

'''


filedir = '/disk/lhcb_data2/RLcMuonic2016/MC/'
inputFile = filedir+'Lb_Lcmunu_MagUp.root'

def main(inputFile):
    F=r.TFile(inputFile)
    T=F.Get("tupleout/DecayTree")

    df0 = r.RDataFrame(T)

    r.gInterpreter.Declare(func_code)
    df1 = df0.Define("isGECPassed" ,  "GECDecision(nVeloClusters, nITClusters,nOTClusters)")
    df2 = df1.Define('isTrackPassed_K',
                     'TrackInputDecision(Lb_PT_DAU_1, Lb_P_DAU_1, Lb_TRACK_CHI2_DAU_1, Lb_TRACK_GHOSTPROB_DAU_1,Lb_TRACK_NDOF_DAU_1)')
    df2 = df2.Define('isTrackPassed_p',
                     'TrackInputDecision(Lb_PT_DAU_2, Lb_P_DAU_2, Lb_TRACK_CHI2_DAU_2, Lb_TRACK_GHOSTPROB_DAU_2,Lb_TRACK_NDOF_DAU_2)')
    df2 = df2.Define('isTrackPassed_pi',
                     'TrackInputDecision(Lb_PT_DAU_3, Lb_P_DAU_3, Lb_TRACK_CHI2_DAU_3, Lb_TRACK_GHOSTPROB_DAU_3,Lb_TRACK_NDOF_DAU_3)')

    df3 = df2.Define('isAdditional_K','TrackAdditionalDecision(K_TRACK_nTTHits)')
    df3 = df3.Define('isAdditional_p','TrackAdditionalDecision(p_TRACK_nTTHits)')
    df3 = df3.Define('isAdditional_pi','TrackAdditionalDecision(pi_TRACK_nTTHits)')

    df4 = df3.Define('isPassed_K','TrackMVADecision(Lb_PT_DAU_1, Lb_TRACK_CHI2_DAU_1, Lb_TRACK_GHOSTPROB_DAU_1,\
                      Lb_IPCHI2_OWNPV_DAU_1,Lb_TRACK_NDOF_DAU_1)')
    df4 = df4.Define('isPassed_p','TrackMVADecision(Lb_PT_DAU_2, Lb_TRACK_CHI2_DAU_2, Lb_TRACK_GHOSTPROB_DAU_2,\
                      Lb_IPCHI2_OWNPV_DAU_2,Lb_TRACK_NDOF_DAU_2)')
    df4 = df4.Define('isPassed_pi','TrackMVADecision(Lb_PT_DAU_3, Lb_TRACK_CHI2_DAU_3, Lb_TRACK_GHOSTPROB_DAU_3,\
                      Lb_IPCHI2_OWNPV_DAU_3,Lb_TRACK_NDOF_DAU_3)')

    df5 = df4.Define('isTrackReco_K','TrackReconstructedDecision()')
    df5 = df5.Define('isTrackReco_p','TrackReconstructedDecision()')
    df5 = df5.Define('isTrackReco_pi','TrackReconstructedDecision()')

    df6 = df5.Define('HLT1TrackMVA_TOS_K','HLT1TrackMVA_TOS(isGECPassed, isTrackPassed_K, isAdditional_K, isPassed_K)')
    df6 = df6.Define('HLT1TrackMVA_TOS_p','HLT1TrackMVA_TOS(isGECPassed, isTrackPassed_p, isAdditional_p, isPassed_p)')
    df6 = df6.Define('HLT1TrackMVA_TOS_pi','HLT1TrackMVA_TOS(isGECPassed, isTrackPassed_pi, isAdditional_pi, isPassed_pi)')

    df7 = df6.Define('HLT1TrackMVA_TOS_Corr_K','HLT1TrackMVA_TOS_EffCorrected(isGECPassed,isTrackPassed_K, isAdditional_K, isPassed_K, isTrackReco_K)')
    df7 = df7.Define('HLT1TrackMVA_TOS_Corr_p','HLT1TrackMVA_TOS_EffCorrected(isGECPassed, isTrackPassed_p, isAdditional_p, isPassed_p, isTrackReco_p)')
    df7 = df7.Define('HLT1TrackMVA_TOS_Corr_pi','HLT1TrackMVA_TOS_EffCorrected(isGECPassed, isTrackPassed_pi, isAdditional_pi, isPassed_pi, isTrackReco_pi)')

    df8 = df7.Define('Lc_HLT1TrackMVA_Emu_TOS','EmulationHLT1(HLT1TrackMVA_TOS_K,HLT1TrackMVA_TOS_p, HLT1TrackMVA_TOS_pi)')
    df8 = df8.Define('Lc_HLT1TrackMVA_Emu_EffCorrected_TOS','EmulationHLT1(HLT1TrackMVA_TOS_Corr_K,HLT1TrackMVA_TOS_Corr_p, HLT1TrackMVA_TOS_Corr_pi)')

    column_name_vector = r.std.vector('string')()
    '''
    column_name_vector.push_back("isGECPassed")
    column_name_vector.push_back("isTrackPassed_K")
    column_name_vector.push_back("isTrackPassed_p")
    column_name_vector.push_back("isTrackPassed_pi")
    column_name_vector.push_back("isAdditional_K")
    column_name_vector.push_back("isAdditional_p")
    column_name_vector.push_back("isAdditional_pi")
    column_name_vector.push_back("isPassed_K")
    column_name_vector.push_back("isPassed_p")
    column_name_vector.push_back("isPassed_pi")
    column_name_vector.push_back("isTrackReco_K")
    column_name_vector.push_back("isTrackReco_p")
    column_name_vector.push_back("isTrackReco_pi")
    column_name_vector.push_back("HLT1TrackMVA_TOS_K")
    column_name_vector.push_back("HLT1TrackMVA_TOS_p")
    column_name_vector.push_back("HLT1TrackMVA_TOS_pi")
    column_name_vector.push_back("HLT1TrackMVA_TOS_Corr_K")
    column_name_vector.push_back("HLT1TrackMVA_TOS_Corr_p")
    column_name_vector.push_back("HLT1TrackMVA_TOS_Corr_pi")
    '''
    column_name_vector.push_back("Lc_HLT1TrackMVA_Emu_TOS")
    column_name_vector.push_back("Lc_HLT1TrackMVA_Emu_EffCorrected_TOS")

    print("Writing the new tree ...")
    new_tree = r.std.vector('string')()
    new_tree.push_back("")
    df8.Snapshot("DecayTree",inputFile[0:-5]+"_wHLT1OneTrackEmulation_Df.root",column_name_vector)
    print("Tree written.")
