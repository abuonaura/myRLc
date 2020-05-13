import ROOT as r
import math as m
import sys, os
import math 
from math import log
import random



def TwoTrackInputDecision(variables):
        if (variables["PT"] <= 600.): return False
        if (variables["P"]  <= 5000.): return False
        if (variables["TRCHI2"] >= 2.5 ): return False
        if (variables["GhostProb"] >= 0.2): return False 
        if (variables["IP"]) <=4.0: return False
        return True


def GetTwoTrackVariables(ev, combination):
        variables={}

        variables["VCHI2"]                      = getattr(ev, "Lb_VERTEX_CHI2_COMB_"+combination)
        variables["ETA"]                        = getattr(ev, "Lb_ETA_COMB_"+combination)
        variables["BPVMCORR"]           = getattr(ev, "Lb_MCORR_OWNPV_COMB_"+combination) 
        variables["BPVDIRA"]        = getattr(ev, "Lb_DIRA_OWNPV_COMB_"+combination)
        variables["MVA"]            = getattr(ev, "Lb_HLt1TwoTrackMVAEmulations_"+combination)
        variables["DOCA"]           = getattr(ev, "Lb_DOCA_COMB_"+combination)
        variables["VDCHI2"]         = getattr(ev, "Lb_VDCHI2_OWNPV_COMB_"+combination)
        
        first_daughter_index  = combination[0]
        second_daughter_index = combination[2]

        if (first_daughter_index == "1"):
                first_daughter_name = "K"
        elif (first_daughter_index == "4"):
                first_daughter_name = "mu"
                
        elif (first_daughter_index == "2"):
                first_daughter_name = "p"
        elif (first_daughter_index == "3"):
                first_daughter_name = "pi"


        if (second_daughter_index == "1"):
                second_daughter_name = "K"
        elif (second_daughter_index == "4"):
                second_daughter_name = "mu"

        elif (second_daughter_index == "2"):
                second_daughter_name = "p"

        elif (second_daughter_index == "3"):
                second_daughter_name = "pi"
                
                

        PX = getattr(ev, first_daughter_name+"_PX") + getattr(ev, second_daughter_name+"_PX")
        PY = getattr(ev, first_daughter_name+"_PY") + getattr(ev, second_daughter_name+"_PY")
        PZ = getattr(ev, first_daughter_name+"_PZ") + getattr(ev, second_daughter_name+"_PZ")

        ASUMPT = math.sqrt(PX*PX + PY*PY)
        ASUMP  = math.sqrt(PX*PX + PY*PY + PZ*PZ)

        variables["SUMPT"] = ASUMPT
        variables["SUMP"]  = ASUMP
        
        variables["PT_COMB"] = getattr(ev, "Lb_PT_COMB_"+combination)
        variables["P_COMB"]  = getattr(ev, "Lb_P_COMB_"+combination) 
        variables["IP"]           = getattr(ev, "Lb_IPCHI2_OWNPV_COMB_"+combination) 

        return variables




def Inefficiency_Decision(Correction):
        #TODO: Check how random is this generator
        rand = random.uniform(0.,1.)
        if rand > Correction: return False
        
        return True     


def TrackReconstructedDecision(variables):
        if(variables["IsReconstructed"] == 0 ): return False
        return True


def TrackMVADecision(variables):
        if ((variables["PT"]<500) or (variables["TRCHI2"]>2.5)):
                return False
        if ((variables["GhostProb"]>=0.2)):
                return False
        if (variables["PT"]>25000):
                if(variables["IP"]>7.4): 
                        return True
                else:
                        return False
        if (variables["IP"]>0) and (variables["PT"]!=500):
#               if (log(variables["IP"]) > 1.0/ pow(variables["PT"]/1000. - 1.0,2) + 1.1/25000.*(25000. - variables["PT"]) + log(7.4) ):
                if (log(variables["IP"]) > 1.0/ pow(variables["PT"]/1000. - 1.0,2) + 1.1/25000.*(25000. - variables["PT"]) + log(7.4) ):
                        return True
        return False


def TwoTrackMVADecision(variables):
        if (variables["VCHI2"] >= 10.): return False
        if (variables["VCHI2"] <= 0. ): return False
        if (variables["ETA"]    <=2. )  : return False
        if (variables["ETA"]   >=5. )  : return False
        if (variables["BPVMCORR"] <= 1000.): return False
        if (variables["BPVMCORR"] >= 1000000000.0): return False
        if (variables["BPVMCORR"] <=0.) : return False
        if (variables["BPVDIRA"]  <= 0.): return False
        if (variables["MVA"]      <= 0.95): return False
        if (variables["SUMPT"]    <= 2000): return False
        if (variables["DOCA"]     >= 10.):  return False
        if (variables["DOCA"]     <= 0.) : return False

        # Remove dummmy values from the input to the MVA
        if(variables["VDCHI2"] <= 0): return False
        if(variables["SUMPT"] <= 0): return False

        if(variables["PT_COMB"] <= 500):  return False
        if(variables["P_COMB"] <= 0.): return False
        if(variables["IP"] <= 4): return False
        if(variables["IP"] < 0): return False 

        return True



def TrackAdditionalDecision(variables):
        if (variables["nTTHits"] < 3):   return False
        return True

def TrackInputDecision(variables):
        if (variables["PT"] <= 500):  return False
        if (variables["P"]  <= 5000): return False
        if (variables["GhostProb"] > 999.): return False
        if (variables["TRCHI2"] >= 4.0): return False
        return True

def GetVariables(head, ev):
        variables={}
        variables["PT"]                  = getattr(ev, "Lb"+"_PT_"+head)
        variables["P"]   = getattr(ev, "Lb"+"_P_"+head)
        variables["IP"]  = getattr(ev, "Lb"+"_IPCHI2_OWNPV_"+head)
        variables["GhostProb"] = getattr(ev, "Lb"+"_TRACK_GHOSTPROB_"+head)
        variables["TRCHI2"]   = float(getattr(ev, "Lb"+"_TRACK_CHI2_"+head)) / getattr(ev, "Lb"+"_TRACK_NDOF_"+head)
        return variables


def AdditionalVariables(ev, head, Correction):
        variables={}
        if head == "DAU_1":
                variables["nTTHits"] = getattr(ev, "K_TRACK_nTTHits")
                variables["History"] = getattr(ev, "K_TRACK_History")
        if head == "DAU_4":
                variables["nTTHits"] = getattr(ev, "mu_TRACK_nTTHits")
                variables["History"] = getattr(ev, "mu_TRACK_History")
        if head == "DAU_2":
                variables["nTTHits"] = getattr(ev, "p_TRACK_nTTHits")
                variables["History"] = getattr(ev, "p_TRACK_History")
        if head == "DAU_3":
                variables["nTTHits"] = getattr(ev, "pi_TRACK_nTTHits")
                variables["History"] = getattr(ev, "pi_TRACK_History")
        variables["IsReconstructed"] =  Inefficiency_Decision(Correction)

        return variables



def GetGECVariables(ev):

        variables={}
        variables["nVeloClusters"] = getattr(ev, "nVeloClusters")
        variables["nITClusters"]   = getattr(ev, "nITClusters")
        variables["nOTClusters"]   = getattr(ev, "nOTClusters")
        variables["nPVs"]          = getattr(ev, "nPVs")

        return variables


def GECDecision(variables):
        if (variables["nVeloClusters"]>=  6000): return False
        if (variables["nOTClusters"]  >= 15000): return False
        if (variables["nITClusters"]  >=  3000): return False
        if (variables["nITClusters"]  <=    50): return False
        if (variables["nOTClusters"]  <=    50): return False
        if (variables["nVeloClusters"] <=   50): return False
        return True

def main1TrackMVA(_name = "Lb_Lcmunu_MagUp.root"):

    #---Get the Track Reconstruction Efficiency correction factor
    Correction = 1./1.035 #This is a number I get from LHCb-PUB-2015-024

    #---Set the random number generator, but seed it for reproducibility
    random.seed(0)

    #_name = "tupleoutMC_trackeronly_v19"
    fin = r.TFile.Open(_name, "READ")
    tin = fin.Get("tupleout/DecayTree")

    for br in ["nVeloClusters", "nITClusters", "nOTClusters", "nPVs"]:
            tin.SetBranchStatus(br, 1)

    for br in ["Lb_TRACK_NDOF_DAU*", "Lb_PT_DAU*", "Lb_P_DAU*", "Lb_IPCHI2_OWNPV_DAU*", "Lb_TRACK_GHOSTPROB_DAU*", "Lb_TRACK_CHI2_DAU*"]:
            tin.SetBranchStatus(br, 1)

    for br in ["*COMB*", "Lb_HLt1TwoTrackMVAEmulations*"]:
            tin.SetBranchStatus(br,1)


    for head in ["mu", "K", "pi", "p"]:
            for br in ["PT", "IPCHI2_OWNPV","PX", "PY", "PZ", "PT", "P"]:
                    tin.SetBranchStatus(head+"_"+br,1)


                    
    fout = r.TFile.Open(_name[:-5]+'_wHLT1OneTrackEmulation.root','RECREATE')
    tout = r.TTree("DecayTree", "DecayTree")

    r.gROOT.ProcessLine("struct MyStruct0{Bool_t abool;};")
    r.gROOT.ProcessLine("struct MyStruct{Int_t aint;};")
    r.gROOT.ProcessLine("struct MyStruct2{Float_t afloat;};")
    from ROOT import MyStruct0
    from ROOT import MyStruct
    from ROOT import MyStruct2

    
    isGECPassed= MyStruct0()
    isTrackPassed_K= MyStruct0()
    isTrackPassed_p= MyStruct0()
    isTrackPassed_pi= MyStruct0()
    isAdditional_K= MyStruct0()
    isAdditional_p= MyStruct0()
    isAdditional_pi= MyStruct0()
    isPassed_K= MyStruct0()
    isPassed_p= MyStruct0()
    isPassed_pi= MyStruct0()
    isTrackReco_K= MyStruct0()
    isTrackReco_p= MyStruct0()
    isTrackReco_pi= MyStruct0()
    HLT1TrackMVA_TOS_K = MyStruct0()
    HLT1TrackMVA_TOS_p = MyStruct0()
    HLT1TrackMVA_TOS_pi = MyStruct0()
    HLT1TrackMVA_TOS_Corr_K = MyStruct0()
    HLT1TrackMVA_TOS_Corr_p = MyStruct0()
    HLT1TrackMVA_TOS_Corr_pi = MyStruct0()


    Lc_HLT1TrackMVA_Emu_TOS                             = MyStruct()
    Lc_HLT1TrackMVA_Emu_EffCorrected_TOS = MyStruct() 

    Lc_HLT1TwoTrackMVA_Emu_TOS              = MyStruct()
    Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS = MyStruct()



    mu_Phi  = MyStruct2()
    K_Phi   = MyStruct2()
    pi_Phi = MyStruct2()
    p_Phia = MyStruct2()

    mu_Theta  = MyStruct2()
    K_Theta   = MyStruct2()
    pi_Theta = MyStruct2()
    p_Theta  = MyStruct2()

    isGECPassed_branch = tout.Branch("isGECPassed",r.AddressOf(isGECPassed,'abool'), "isGECPassed/B")
    
    isTrackPassed_branch_K = tout.Branch("isTrackPassed_K",r.AddressOf(isTrackPassed_K,'abool'), "isTrackPassed_K/B")
    isTrackPassed_branch_p = tout.Branch("isTrackPassed_p",r.AddressOf(isTrackPassed_p,'abool'), "isTrackPassed_p/B")
    isTrackPassed_branch_pi = tout.Branch("isTrackPassed_pi",r.AddressOf(isTrackPassed_pi,'abool'), "isTrackPassed_pi/B")

    isAdditional_branch_K = tout.Branch("isAdditional_K",r.AddressOf(isAdditional_K,'abool'), "isAdditional_K/B")
    isAdditional_branch_p = tout.Branch("isAdditional_p",r.AddressOf(isAdditional_p,'abool'), "isAdditional_p/B")
    isAdditional_branch_pi = tout.Branch("isAdditional_pi",r.AddressOf(isAdditional_pi,'abool'), "isAdditional_pi/B")
    
    isPassed_branch_K = tout.Branch("isPassed_K",r.AddressOf(isPassed_K,'abool'), "isPassed_K/B")
    isPassed_branch_p = tout.Branch("isPassed_p",r.AddressOf(isPassed_p,'abool'), "isPassed_p/B")
    isPassed_branch_pi = tout.Branch("isPassed_pi",r.AddressOf(isPassed_pi,'abool'), "isPassed_pi/B")
    
    isTrackReco_branch_K = tout.Branch("isTrackReco_K",r.AddressOf(isTrackReco_K,'abool'), "isTrackReco_K/B")
    isTrackReco_branch_p = tout.Branch("isTrackReco_p",r.AddressOf(isTrackReco_p,'abool'), "isTrackReco_p/B")
    isTrackReco_branch_pi = tout.Branch("isTrackReco_pi",r.AddressOf(isTrackReco_pi,'abool'), "isTrackReco_pi/B")

    HLT1TrackMVA_Emu_TOS_branch_K = tout.Branch("HLT1TrackMVA_TOS_K",r.AddressOf(HLT1TrackMVA_TOS_K,'abool'), "HLT1TrackMVA_TOS_K/B")
    HLT1TrackMVA_Emu_TOS_branch_p = tout.Branch("HLT1TrackMVA_TOS_p",r.AddressOf(HLT1TrackMVA_TOS_p,'abool'), "HLT1TrackMVA_TOS_p/B")
    HLT1TrackMVA_Emu_TOS_branch_pi = tout.Branch("HLT1TrackMVA_TOS_pi",r.AddressOf(HLT1TrackMVA_TOS_pi,'abool'), "HLT1TrackMVA_TOS_pi/B")

    HLT1TrackMVA_Emu_TOS_Corr_branch_K = tout.Branch("HLT1TrackMVA_TOS_Corr_K",r.AddressOf(HLT1TrackMVA_TOS_Corr_K,'abool'), "HLT1TrackMVA_TOS_Corr_K/B")
    HLT1TrackMVA_Emu_TOS_Corr_branch_p = tout.Branch("HLT1TrackMVA_TOS_Corr_p",r.AddressOf(HLT1TrackMVA_TOS_Corr_p,'abool'), "HLT1TrackMVA_TOS_Corr_p/B")
    HLT1TrackMVA_Emu_TOS_Corr_branch_pi = tout.Branch("HLT1TrackMVA_TOS_Corr_pi",r.AddressOf(HLT1TrackMVA_TOS_Corr_pi,'abool'), "HLT1TrackMVA_TOS_Corr_pi/B")
    
    Lc_HLT1TrackMVA_Emu_TOS_branch = tout.Branch("Lc_HLT1TrackMVA_Emu_TOS",r.AddressOf(Lc_HLT1TrackMVA_Emu_TOS,'aint'), "Lc_HLT1TrackMVA_Emu_TOS/I")
    Lc_HLT1TrackMVA_Emu_EffCorrected_TOS_branch = tout.Branch("Lc_HLT1TrackMVA_Emu_EffCorrected_TOS", r.AddressOf(Lc_HLT1TrackMVA_Emu_EffCorrected_TOS, 'aint'), "Lc_HLT1TrackMVA_Emu_EffCorrected_TOS/I") 

    Lc_HLT1TwoTrackMVA_Emu_TOS_branch = tout.Branch("Lc_HLT1TwoTrackMVA_Emu_TOS",r.AddressOf(Lc_HLT1TwoTrackMVA_Emu_TOS,'aint'),"Lc_HLT1TwoTrackMVA_Emu_TOS/I")
    Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS_branch = tout.Branch("Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS", r.AddressOf(Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS, 'aint'), "Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS/I")



    mu_Phi_branch  = tout.Branch("mu_Phi",  r.AddressOf(mu_Phi,  "afloat"), "mu_Phi/F")
    K_Phi_branch   = tout.Branch("K_Phi",   r.AddressOf(K_Phi,   "afloat"), "K_Phi/F")
    pi_Phi_branch = tout.Branch("pi_Phi", r.AddressOf(pi_Phi, "afloat"), "pi_Phi/F")
    p_Phi_branch = tout.Branch("p_Phi", r.AddressOf(p_Phia, "afloat"), "p_Phi/F")

    mu_Theta_branch  = tout.Branch("mu_Theta",  r.AddressOf(mu_Theta,  "afloat"), "mu_Theta/F")
    K_Theta_branch   = tout.Branch("K_Theta",   r.AddressOf(K_Theta,   "afloat"), "K_Theta/F")
    pi_Theta_branch = tout.Branch("pi_Theta", r.AddressOf(pi_Theta, "afloat"), "pi_Theta/F")
    p_Theta_branch = tout.Branch("p_Theta", r.AddressOf(p_Theta, "afloat"), "p_Theta/F")


    #-------------------------------------------------------------------------------------------------------------

    NEvents = tin.GetEntries()

    i=0
    print ("Processing events ...")
    for ev in tin:
            if (i%1000 == 0):
                    sys.stdout.write("\rProcessed "+str(round(float(i)/NEvents*100, 2))+"%")
                    sys.stdout.flush()
            #if (i == 15000):
            #    break
                
            Lc_HLT1TrackMVA_Emu_TOS.aint                                = 0 
            Lc_HLT1TrackMVA_Emu_EffCorrected_TOS.aint   = 0

            Lc_HLT1TwoTrackMVA_Emu_TOS.aint             = 0
            Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS.aint = 0

            mu_Phi.afloat  = 0
            K_Phi.afloat   = 0
            pi_Phi.afloat = 0
            p_Phia.afloat = 0

            mu_Theta.afloat  = 0
            K_Theta.afloat       = 0
            pi_Theta.afloat = 0
            p_Theta.afloat = 0
            
            #------------------------------------------------------------------------------
            
            
            #---Global event decisions
            GECvariables = GetGECVariables(ev)
            IsGECPassed  = GECDecision(GECvariables)
            isGECPassed.abool = IsGECPassed
            
                    #---Save all needed informations once for all
            Variables_vec                       = {}
            AdditionalVariables_vec = {}
            TwoTrackVariables_vec   = {}


            for head in ["DAU_1", "DAU_2", "DAU_3", "DAU_4"]:
                    Variables_vec[head]                   = GetVariables(head, ev)
                    AdditionalVariables_vec[head] = AdditionalVariables(ev, head, Correction)
        
            for head in ["1_2", "1_3", "1_4", "2_3", "2_4", "3_4"]:
                    TwoTrackVariables_vec[head] = GetTwoTrackVariables(ev, head)



            #---Single particle decisions
            for head in ["DAU_1", "DAU_2", "DAU_3"]:
                    IsTrackPassed                 = TrackInputDecision(Variables_vec[head])
                    IsAdditional              = TrackAdditionalDecision(AdditionalVariables_vec[head])  
                    IsPassed                      = TrackMVADecision(Variables_vec[head])
                    IsTrackReconstructed  = TrackReconstructedDecision(AdditionalVariables_vec[head])
                    if head=="DAU_1":
                        isTrackPassed_K.abool = IsTrackPassed 
                        isAdditional_K.abool = IsAdditional
                        isPassed_K.abool = IsPassed 
                        isTrackReco_K.abool = IsTrackReconstructed
                    if head=="DAU_2":
                        isTrackPassed_p.abool = IsTrackPassed 
                        isAdditional_p.abool = IsAdditional
                        isPassed_p.abool = IsPassed 
                        isTrackReco_p.abool = IsTrackReconstructed
                    if head=="DAU_3":
                        isTrackPassed_pi.abool = IsTrackPassed 
                        isAdditional_pi.abool = IsAdditional
                        isPassed_pi.abool = IsPassed 
                        isTrackReco_pi.abool = IsTrackReconstructed

                    
                    if (IsTrackReconstructed and IsPassed and IsGECPassed and IsTrackPassed and IsAdditional):
                            Lc_HLT1TrackMVA_Emu_EffCorrected_TOS.aint = 1
                    if (IsPassed and IsGECPassed and IsTrackPassed and IsAdditional):
                            Lc_HLT1TrackMVA_Emu_TOS.aint = 1

            '''


            #---Two particle decisions
            for k in range(1,3):
                    j = k+1
                    while (j<=3):

                            combination  = str(k)+"_"+str(j)

                            variables_k           = Variables_vec["DAU_"+str(k)]
                            variables_j           = Variables_vec["DAU_"+str(j)]

                            AdditionalVariables_k = AdditionalVariables_vec["DAU_"+str(k)]
                            AdditionalVariables_j = AdditionalVariables_vec["DAU_"+str(j)]

                            IsTrackPassed_k       = TrackInputDecision(variables_k)
                            IsTrackPassed_j       = TrackInputDecision(variables_j)

                            IsTrackPassed  = IsTrackPassed_k and IsTrackPassed_j
                    
                            IsAdditional_k        = TrackAdditionalDecision(AdditionalVariables_k)
                            IsAdditional_j        = TrackAdditionalDecision(AdditionalVariables_j)

                            IsAdditional   = IsAdditional_k and IsAdditional_j

                            IsInputPassed_k = TwoTrackInputDecision(variables_k)
                            IsInputPassed_j = TwoTrackInputDecision(variables_j)

                            IsInputPassed                 = IsInputPassed_k and IsInputPassed_j

                            variables_comb                = TwoTrackVariables_vec[combination]

                            IsPassed              = TwoTrackMVADecision(variables_comb)

                            IsTrackReconstructed_k  = TrackReconstructedDecision(AdditionalVariables_k)
                            IsTrackReconstructed_j  = TrackReconstructedDecision(AdditionalVariables_j)

                            IsTrackReconstructed = IsTrackReconstructed_k and IsTrackReconstructed_j

                            if (IsPassed                    and IsInputPassed              and IsGECPassed and IsTrackPassed and IsAdditional):
                                    Lc_HLT1TwoTrackMVA_Emu_TOS.aint = 1

                            if (IsTrackReconstructed and IsPassed                and IsInputPassed  and IsGECPassed and IsTrackPassed and IsAdditional):
                                    Lc_HLT1TwoTrackMVA_Emu_EffCorrected_TOS.aint = 1

                            j = j+1

            #-----Fill the additional variables
            mu_Phi.afloat  = math.atan(getattr(tin, "mu_PY")  / getattr(tin, "mu_PX" ))
            K_Phi.afloat   = math.atan(getattr(tin, "K_PY")   / getattr(tin, "K_PX"  ))
            p_Phia.afloat = math.atan(getattr(tin, "p_PY") / getattr(tin, "p_PX"))
            pi_Phi.afloat = math.atan(getattr(tin, "pi_PY") / getattr(tin, "pi_PX"))

            mu_Theta.afloat  = math.acos(getattr(tin, "mu_PZ")  / getattr(tin, "mu_P" ))
            K_Theta.afloat   = math.acos(getattr(tin, "K_PZ")   / getattr(tin, "K_P"  ))
            p_Theta.afloat = math.acos(getattr(tin, "p_PZ") / getattr(tin, "p_P"))
            pi_Theta.afloat = math.acos(getattr(tin, "pi_PZ") / getattr(tin, "pi_P"))
            '''

            i = i+1
            tout.Fill()
    print ("\nAll events processed.")
    tout.Write()
    fout.Close()
    fin.Close()


        
        






