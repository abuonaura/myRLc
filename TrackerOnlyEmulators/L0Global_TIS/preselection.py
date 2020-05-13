from math import log
import ROOT

def pseudorapidity(p_,pz_):
   return 0.5*log((p_+pz_)/(p_-pz_))

ROOT.gROOT.ProcessLine(\
	"struct MyStruct{\
	Float_t afloat;\
	};")
from ROOT import MyStruct

# Truth matching
cut_string = ""

############
# J/psi cuts

cut_string += "mup_PIDmu>0"
cut_string += " && mum_PIDmu>0"
cut_string += " && mup_PT>500"
cut_string += " && mum_PT>500"
#cut_string += " && Jpsi_TAU>0.0002"
cut_string += " && abs(Jpsi_MM-3096.9)<80"
cut_string += " && Jpsi_ENDVERTEX_CHI2<16"

############
# K cuts

cut_string += " && Kplus_IPCHI2_OWNPV>25"
cut_string += " && Kplus_PIDK>0"
cut_string += " && Kplus_TRACK_GhostProb<0.35"

############
# Bp cuts

cut_string += " && (Bplus_ENDVERTEX_CHI2/Bplus_ENDVERTEX_NDOF)<6"
cut_string += " && Bplus_M>(5280-50)"
cut_string += " && Bplus_M<(5280+50)"

############
# Trigger

#cut_string += " && Bplus_L0MuonDecision_TOS" # Added at a later stage
#cut_string += " && Bplus_Hlt1TrackMuonDecision_TOS" # Added at a later stage
#cut_string += " && Bplus_Hlt2DiMuonDetachedHeavyDecision_TOS" # Added at a later stage

br_list = []

br_list.append("Bplus_L0Global*")
br_list.append("Bplus_L0MuonDecision_TOS")
br_list.append("Bplus_L0DiMuonDecision_TOS")
br_list.append("Bplus_Hlt1TrackMuonDecision_TOS")
br_list.append("Bplus_Hlt2DiMuonDetachedHeavyDecision_TOS")
br_list.append("mup_PIDmu")
br_list.append("mum_PIDmu")
br_list.append("mup_PT")
br_list.append("mum_PT")
br_list.append("Jpsi_MM")
br_list.append("Jpsi_ENDVERTEX_CHI2")
br_list.append("Kplus_IPCHI2_OWNPV")
br_list.append("Kplus_PIDK")
br_list.append("Kplus_TRACK_GhostProb")
br_list.append("Bplus_ENDVERTEX_CHI2")
br_list.append("Bplus_ENDVERTEX_NDOF")
br_list.append("Bplus_M")
br_list.append("nSPDHits")
br_list.append("Bplus_P")
br_list.append("Bplus_PT")
br_list.append("Bplus_PZ")
br_list.append("Jpsi_P")
br_list.append("Jpsi_PT")
br_list.append("Jpsi_PZ")
br_list.append("Kplus_P")
br_list.append("Kplus_PT")
br_list.append("Kplus_PZ")
br_list.append("*L0*")
br_list.append("*Hlt*")
br_list.append("nTracks")
