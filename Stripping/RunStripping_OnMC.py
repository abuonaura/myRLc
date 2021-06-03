import  GaudiKernel.SystemOfUnits as Units
from Gaudi.Configuration import *
from Configurables import DaVinci
from PhysSelPython.Wrappers import *
import copy
import pdb
import re

##########################################################################
#####                    DA VINCI SETTINGS
#########################################################################

DaVinci().HistogramFile = 'DV_stripping_histos.root'
#DaVinci().HistogramFile = 'DV_stripping_histos_LbLcmunu_fullsim.root'
#DaVinci().HistogramFile = 'DV_stripping_histos_LbLcmunu_trackeronly.root'
DaVinci().EvtMax = -1                         
#DaVinci().EvtMax = 5000
DaVinci().PrintFreq = 5000
DaVinci().DataType  = "2016"
DaVinci().Simulation = True
DaVinci().Lumi =  not DaVinci().Simulation
DaVinci().ProductionType = "Stripping"
DaVinci().InputType = "DST"
#DaVinci().TupleFile = "LbLcmunu_fullsim.root"
#DaVinci().TupleFile = "/disk/lhcb_data2/buonaura/LbLcmunu_MagUp_trackeronly_5000.root"
#DaVinci().TupleFile = "LbLcDs_fullsim.root"

DaVinci().TupleFile = "tupleoutMC.root"
#DaVinci().TupleFile = "tupleoutMC_trackeronly.root"

# change the column size of timing table                                                                             
from Configurables import TimingAuditor, SequencerTimerTool
TimingAuditor().addTool(SequencerTimerTool,name="TIMER")
TimingAuditor().TIMER.NameSize = 60

if DaVinci().EvtMax !=-1:
	#DaVinci().Input=['PFN:/home/hep/buonaura/Analysis/RLc/Datasets/MC/00090120_00000030_1.lctaunu.safestrip.dst']
#------ Lb_Lc2593munu
	#DaVinci().Input=['PFN:/home/hep/buonaura/Analysis/RLc/Datasets/MC/00067107_00000006_1.lctaunu.safestrip.dst']
#------ Lb_Lcmunu full sim
#	DaVinci().Input=['PFN:/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/DSTfiles/Lb_Lcmunu/FullSim/00066927_00000004_1.lctaunu.safestrip.dst']
#------ Lb_Lcmunu MagUp trackeronly
	DaVinci().Input=['PFN:/disk/lhcb_data2/buonaura/00089844_00000124_1.lctaunu.safestrip.dst']
#------ Lb_LcDs full sim
   #DaVinci().Input=['PFN:/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/DSTfiles/Lb_LcDs/FullSim/00067063_00000002_1.lctaunu.safestrip.dst']



###########################################################################
#####              IMPORT STRIPPING LINE
###########################################################################
from StrippingConf.StrippingStream import StrippingStream
from StrippingConf.Configuration import StrippingConf
from StrippingSelections.StrippingSL import StrippingB2DMuForTauMu
#import actual selections
from StrippingSettings.Stripping28.LineConfigDictionaries_Semileptonic import B2DMuForTauMu

#**************************************************************************
# MODIFY STRIPPING CUTS for PIDs. 
#---> PID is not applied to MC because badly modeled.
#---> Corrections with data are applied to MC PID variables
#---> Substitute the tauonic line PID cuts.
B2DMuForTauMu['CONFIG']['PIDmu']         = -99999
B2DMuForTauMu['CONFIG']['KaonPIDK']      = -99999
B2DMuForTauMu['CONFIG']['ProtonPIDp']      = -99999
B2DMuForTauMu['CONFIG']['PionPIDKTight'] =  99999
B2DMuForTauMu['CONFIG']['Hlt2Line'] =  ""

MyStream = StrippingStream("B2DMuNuXANDTAU")
#create a line builder instance of type 'B2DMuForTauMu' and configure it with the 'CONFIG' dictionary.
confB2DMuForTauMu = StrippingB2DMuForTauMu.B2DMuForTauMuconf("B2DMuForTauMu", B2DMuForTauMu['CONFIG'])

#Replace Mu/Pi/K/p Selections with NoPID particles
replaceDictFilterDesktop = { 'MuforB2DMuForTauMu'  : 'Phys/StdAllNoPIDsMuons/Particles'
		             ,'PiforB2DMuForTauMu' : 'Phys/StdAllNoPIDsPions/Particles'
			     ,'KforB2DMuForTauMu'  : 'Phys/StdAllNoPIDsKaons/Particles' 
			     ,'PforB2DMuForTauMu'  : 'Phys/StdAllNoPIDsProtons/Particles' }
replaceDictVoidFilter = { 'SelFilterPhys_StdAllLooseMuons_Particles' : "\n 0<CONTAINS('Phys/StdAllNoPIDsMuons/Particles',True)\n "
		          ,'SelFilterPhys_StdLooseKaons_Particles'   : "\n 0<CONTAINS('Phys/StdAllNoPIDsKaons/Particles',True)\n "
			  ,'SelFilterPhys_StdLoosePions_Particles'   : "\n 0<CONTAINS('Phys/StdAllNoPIDsPions/Particles',True)\n "
			  ,'SelFilterPhys_StdLooseProtons_Particles'   : "\n 0<CONTAINS('Phys/StdAllNoPIDsProtons/Particles',True)\n "
			  }

#Read the line builders and loop over them 
lineList2= confB2DMuForTauMu.lines()
for l in lineList2:
	MyLines = []
	#Select the one for our analysis (b2LcMuXB2DMuForTauMuLine)
	if 'b2LcMuXB2DMuForTauMuLine' in l.name():
		#loop over the members in this line:
		for i in range(len(l._members)):
			seq = l._members[i]
			print (seq.name())
			#Replace particles with NoPID ones:
			if seq.name() in replaceDictFilterDesktop.keys():
				seq.Inputs = [replaceDictFilterDesktop[seq.name()]]
				l._members = l._members[:i] + [seq] + l._members[i+1:]
			if seq.name() in replaceDictVoidFilter.keys():
				seq.Code = replaceDictVoidFilter[seq.name()]
				l._members = l._members[:i] + [seq] + l._members[i+1:]
		MyLines.append(l)
	MyStream.appendLines(MyLines)

locations = []
for lin in MyStream.lines:
	print lin.outputLocation()
	locations.append(lin.outputLocation())

#General configurations for stripping
dstStreams  = [ "Semileptonic"]
stripTESPrefix = 'Strip'
from Configurables import ProcStatusCheck
sc = StrippingConf( Streams = [MyStream] ,
		    MaxCandidates = 20000,
		    AcceptBadEvents = False,
		    BadEventSelection = ProcStatusCheck(),
		    TESPrefix = stripTESPrefix,
		    Verbose = True,
		    DSTStreams = dstStreams)

# So that we do not get all events written out
MyStream.sequence().IgnoreFilterPassed = False

#This deletes stripping objects from previous stripping versions in the DST file                             
from Configurables import EventNodeKiller
eventNodeKiller = EventNodeKiller('Stripkiller')
eventNodeKiller.Nodes = [ '/Event/AllStreams', '/Event/Strip' ]

DaVinci().appendToMainSequence([eventNodeKiller,sc.sequence()])

##############################################################################################
##### HLT1 Emulation
##############################################################################################

relinfoutput = MyStream.outputLocations()[0].replace("Particles", "HLT1Emulation")

def getVars(i, j):
	combination = str(i)+"_"+str(j)
	Vars = {'chi2'   : "RELINFO('"+ relinfoutput +"', 'VERTEX_CHI2_COMB_"+combination+"',0)",
	  'fdchi2' : "RELINFO('"+ relinfoutput +"', 'VDCHI2_MINIPPV_COMB_"+combination+"',0)",
	  'sumpt'  : "RELINFO('"+ relinfoutput +"', 'SUMPT_COMB_"+combination+"',0)",
	  'nlt16'  : "RELINFO('"+ relinfoutput +"', 'NLT_MINIPPV_COMB_"+combination+"',0)"
	 }
	return Vars

def getVars_WrongIP(i, j):
	combination = str(i)+"_"+str(j)
	Vars = {'chi2'   : "RELINFO('"+ relinfoutput +"', 'VERTEX_CHI2_COMB_"+combination+"',0)",
	        'fdchi2' : "RELINFO('"+ relinfoutput +"', 'VDCHI2_OWNPV_COMB_"+combination+"',0)",
		'sumpt'  : "RELINFO('"+ relinfoutput +"', 'SUMPT_COMB_"+combination+"',0)",
		'nlt16'  : "RELINFO('"+ relinfoutput +"', 'NLT_OWNPV_COMB_"+combination+"',0)"
		}
	return Vars

from Configurables import RelInfoHLT1Emulation,AddRelatedInfo

relinfo_HLT1 = AddRelatedInfo("RelInfo_HLT1")
relinfo_HLT1.addTool(RelInfoHLT1Emulation, "RelInfoHLT1Emulation")
relinfo_HLT1.Tool = "RelInfoHLT1Emulation"
relinfo_HLT1.Location = 'HLT1Emulation'
relinfo_HLT1.Inputs = [MyStream.outputLocations()[0]]
hlt1_emu = getattr(relinfo_HLT1, "RelInfoHLT1Emulation")
hlt1_emu.Variables = []
hlt1_emu.nltValue = 16.

DaVinci().appendToMainSequence([relinfo_HLT1])

##############################################################################################                      
##### PRODUCTION OF THE VELO PARTICLE INPUT CONTAINER FOR THE CHARGED ISOLATION TUPLE TOOL   
#############################################################################################                       

#To produce the Velo particles input container (the gain in performance from adding velo tracks is actually pretty negligible...)                                                                                                       

from Configurables import ChargedProtoParticleMaker

name = "Lambdab2Lcmunu"
veloprotos = ChargedProtoParticleMaker(name+"ProtoPMaker")
veloprotos.Inputs = ["Rec/Track/Best"]
veloprotos.Output = "Rec/ProtoP/myProtoPMaker/ProtoParticles"


DaVinci().appendToMainSequence( [ veloprotos ])


from Gaudi.Configuration import *
from Configurables       import ProtoParticleCALOFilter, CombinedParticleMaker,NoPIDsParticleMaker
from CommonParticles.Utils import *


algorithm = NoPIDsParticleMaker('StdNoPIDsVeloPions',  Particle = 'pion',  )
algorithm.Input = "Rec/ProtoP/myProtoPMaker/ProtoParticles"
selector = trackSelector ( algorithm , trackTypes = ['Velo'] )

locations = updateDoD ( algorithm )
DaVinci().appendToMainSequence( [ algorithm ])

##############################################################################################
##### TOOLS AND TRIGGERS AND STRIPPING LINES AND LOKI VARIABLES
##############################################################################################
enablePacking = True


from Configurables import DecayTreeTuple, FitDecayTrees, TupleToolRecoStats, TupleToolTrigger, TupleToolSubMass
from Configurables import TupleToolTISTOS, CondDB, SelDSTWriter, TupleToolL0Calo
from Configurables import TupleToolTrackInfo, TupleToolRICHPid, TupleToolGeometry, TupleToolPid
from Configurables import TupleToolANNPID
from Configurables import TupleToolSLTruth
from Configurables import TupleToolMCTruth
from Configurables import TupleToolPropertime
from Configurables import TupleToolTauMuDiscrVarsLcMassConstraint
from DecayTreeTuple.Configuration import *

tupleB = DecayTreeTuple("tupleout")
tupleB.Inputs = ["Phys/b2LcMuXB2DMuForTauMuLine/Particles"]
tupleB.Decay= "[B- -> ^(Lambda_c+ -> ^p+ ^K- ^pi+) ^mu-]CC"
tupleB.ToolList +=  [
                "TupleToolKinematic"
                , "TupleToolEventInfo"
                , "TupleToolRecoStats"
                , "TupleToolGeometry"
                ,"TupleToolMCBackgroundInfo"#comment out for data      
                ,"TupleToolMCTruth" #comment out for data                                           
                ,"TupleToolSLTruth"
                ,"TupleToolL0Calo" #for Calorimeter info                                             
                ] # Probably need to add many more Tools.                         


tupleB.addBranches ({
        "Lb": "[B- -> (Lambda_c+ -> p+ K- pi+) mu-]CC",
        "Lc" : "[B- -> ^Lambda_c+ mu-]CC",
        "p"  : "[B- -> (Lambda_c+ -> ^p+ K- pi+) mu-]CC",
        "pi" : "[B- -> (Lambda_c+ -> p+ K- ^pi+) mu-]CC",
        "K" : "[B- -> (Lambda_c+ -> p+ ^K- pi+) mu-]CC",
        "mu" : "[B- -> (Lambda_c+ -> p+ K- pi+) ^mu-]CC"
        })

l0_lines = ['L0HadronDecision'
            ,'L0MuonDecision'
            ,'L0ElectronDecision'
            ]

hlt1_lines = ['Hlt1TrackMVADecision'
              ,'Hlt1TwoTrackMVADecision'
             ]

hlt2_lines = ["Hlt2XcMuXForTauB2XcMuDecision",
              "Hlt2XcMuXForTauB2XcFakeMuDecision"
              ]


LoKi_All=tupleB.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_All")
LoKi_All.Variables = {
                'MINIPCHI2' : "MIPCHI2DV(PRIMARY)",
                'MINIP' : "MIPDV(PRIMARY)",
                'IPCHI2_OWNPV' : "BPVIPCHI2()",
                'IP_OWNPV' : "BPVIP()",
                'ghost' : "TRGHP",
                'TRACK_CHI2' : "TRCHI2DOF"
                }

tt_taumudiscrvars_b = tupleB.Lb.addTupleTool("TupleToolTauMuDiscrVarsLcMassConstraint")

LoKi_B=tupleB.Lb.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_B")
LoKi_B.Variables = {
                'DOCAMAX' : "DOCAMAX",
                'TAU' : "BPVLTIME()",
                'DIRA_OWNPV' : "BPVDIRA",
                'FD_CHI2' : "BPVVDCHI2",
                'ENDVERTEX_CHI2' : "VFASPF(VCHI2/VDOF)",
                'Corrected_Mass' : "BPVCORRM"
                }

LoKi_B=tupleB.Lc.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_B")
LoKi_B.Variables = {
                'TAU' : "BPVLTIME()"
                }
'''
tupletistos = tupleB.Lb.addTupleTool("TupleToolTISTOS")
tupletistos.VerboseL0 = True
tupletistos.VerboseHlt1 = False
tupletistos.VerboseHlt2 = True
tupletistos.TriggerList = l0_lines + hlt2_lines

tupletistos1= tupleB.Lc.addTupleTool("TupleToolTISTOS")
tupletistos1.VerboseL0 = True
tupletistos1.VerboseHlt1 = True
tupletistos1.VerboseHlt2 = False
tupletistos1.TriggerList = l0_lines + hlt1_lines
'''

LoKi_muplus=tupleB.mu.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_muplus")
LoKi_muplus.Variables = {
                'PIDmu' : "PIDmu",
                'NNmu' : "PPINFO(PROBNNmu)"
		}

tt_trackinfo = tupleB.addTupleTool("TupleToolTrackInfo")
tt_trackinfo.Verbose = True


DaVinci().appendToMainSequence( [tupleB] )

#*************************************                        
#Isolation TupleTool:  
#*************************************

from Configurables import TupleToolApplyIsolationMC
tupleB.Lb.addTool(TupleToolApplyIsolationMC, name="TupleToolApplyIsolationMC")
tupleB.Lb.TupleToolApplyIsolationMC.WeightsFile="weights.xml"
tupleB.Lb.ToolList+=["TupleToolApplyIsolationMC/TupleToolApplyIsolationMC"]

             
#*************************************                        
#TupleTool for Lc lifetime  
#************************************* 
tupleB.Lc.ToolList += [ "TupleToolPropertime" ]

#************************************* 
#TupleTool to save history of particles
#*************************************           
#tupleB.TupleToolMCTruth.ToolList += ["MCTupleToolHierarchy","MCTupleToolHierarchy"]          
MCTruth = TupleToolMCTruth()
MCTruth.ToolList += ["MCTupleToolHierarchy", "MCTupleToolKinematic"]
tupleB.addTool(MCTruth)




############################################################             
# Hlt1 Emulation                                                                                                   
############################################################                                                 

from MVADictHelpers import addMatrixnetclassifierTuple
branch = tupleB.Lb
#branch = tupleB.Lc
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,2), "HLt1TwoTrackMVAEmulations_1_2", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,3), "HLt1TwoTrackMVAEmulations_1_3", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(1,4), "HLt1TwoTrackMVAEmulations_1_4", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(2,3), "HLt1TwoTrackMVAEmulations_2_3", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(2,4), "HLt1TwoTrackMVAEmulations_2_4", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars(3,4), "HLt1TwoTrackMVAEmulations_3_4", True)

addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(1,2), "HLt1TwoTrackMVAEmulations_WrongIP_1_2", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(1,3), "HLt1TwoTrackMVAEmulations_WrongIP_1_3", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(1,4), "HLt1TwoTrackMVAEmulations_WrongIP_1_4", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(2,3), "HLt1TwoTrackMVAEmulations_WrongIP_2_3", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(2,4), "HLt1TwoTrackMVAEmulations_WrongIP_2_4", True)
addMatrixnetclassifierTuple(branch, "$PARAMFILESROOT/../v8r27p1/data/Hlt1TwoTrackMVA.mx", getVars_WrongIP(3,4), "HLt1TwoTrackMVAEmulations_WrongIP_3_4", True)

# Hlt1 emulation variables                                                                                        
tt_HLT1Emulation = branch.addTupleTool("LoKi::Hybrid::TupleTool/Hlt1TwoTrackMVAEmulation")
tt_HLT1Emulation.Preambulo = []
tt_HLT1Emulation.Variables = {
                "NDAUGHTERS" :       "RELINFO('"+ relinfoutput+ "', 'NDAUGHTERS', 0)",
                "PT_DAU_1"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_1', 0)",
                "PT_DAU_2"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_2', 0)",
                "PT_DAU_3"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_3', 0)",
                "PT_DAU_4"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_4', 0)",
                "PT_DAU_5"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_5', 0)",
                "PT_DAU_6"   :       "RELINFO('"+ relinfoutput+ "', 'PT_DAU_6', 0)",

                "P_DAU_1"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_1', 0)",
                "P_DAU_2"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_2', 0)",
                "P_DAU_3"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_3', 0)",
                "P_DAU_4"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_4', 0)",
                "P_DAU_5"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_5', 0)",
                "P_DAU_6"   :       "RELINFO('"+ relinfoutput+ "', 'P_DAU_6', 0)",

                "IPCHI2_OWNPV_DAU_1" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_1', 0)",
                "IPCHI2_OWNPV_DAU_2" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_2', 0)",
                "IPCHI2_OWNPV_DAU_3" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_3', 0)",
                "IPCHI2_OWNPV_DAU_4" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_4', 0)",
                "IPCHI2_OWNPV_DAU_5" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_5', 0)",
                "IPCHI2_OWNPV_DAU_6" : "RELINFO('"+ relinfoutput+ "', 'IPCHI2_OWNPV_DAU_6', 0)",

                "TRACK_GHOSTPROB_DAU_1" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_1', 0)",
                "TRACK_GHOSTPROB_DAU_2" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_2', 0)",
                "TRACK_GHOSTPROB_DAU_3" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_3', 0)",
                "TRACK_GHOSTPROB_DAU_4" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_4', 0)",
                "TRACK_GHOSTPROB_DAU_5" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_5', 0)",
                "TRACK_GHOSTPROB_DAU_6" : "RELINFO('"+ relinfoutput+ "', 'TRACK_GHOSTPROB_DAU_6', 0)",

                "TRACK_CHI2_DAU_1" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_1', 0)",
                "TRACK_CHI2_DAU_2" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_2', 0)",
                "TRACK_CHI2_DAU_3" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_3', 0)",
                "TRACK_CHI2_DAU_4" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_4', 0)",
                "TRACK_CHI2_DAU_5" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_5', 0)",
                "TRACK_CHI2_DAU_6" : "RELINFO('"+ relinfoutput+ "', 'TRACK_CHI2_DAU_6', 0)",

                "TRACK_NDOF_DAU_1" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_1', 0)",
                "TRACK_NDOF_DAU_2" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_2', 0)",
                "TRACK_NDOF_DAU_3" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_3', 0)",
                "TRACK_NDOF_DAU_4" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_4', 0)",
                "TRACK_NDOF_DAU_5" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_5', 0)",
                "TRACK_NDOF_DAU_6" : "RELINFO('"+ relinfoutput+ "', 'TRACK_NDOF_DAU_6', 0)",

		"VERTEX_CHI2_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_1_2',0)",
                "VERTEX_CHI2_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_1_3',0)",
                "VERTEX_CHI2_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_1_4',0)",
                "VERTEX_CHI2_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_1_5',0)",
                "VERTEX_CHI2_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_1_6',0)",
                "VERTEX_CHI2_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_2_3',0)",
                "VERTEX_CHI2_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_2_4',0)",
                "VERTEX_CHI2_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_2_5',0)",
                "VERTEX_CHI2_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_2_6',0)",
                "VERTEX_CHI2_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_3_4',0)",
                "VERTEX_CHI2_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_3_5',0)",
                "VERTEX_CHI2_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_3_6',0)",
                "VERTEX_CHI2_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_4_5',0)",
                "VERTEX_CHI2_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_4_6',0)",
                "VERTEX_CHI2_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_CHI2_COMB_5_6',0)",

                "VERTEX_NDOF_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_1_2',0)",
                "VERTEX_NDOF_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_1_3',0)",
                "VERTEX_NDOF_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_1_4',0)",
                "VERTEX_NDOF_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_1_5',0)",
                "VERTEX_NDOF_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_1_6',0)",
                "VERTEX_NDOF_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_2_3',0)",
                "VERTEX_NDOF_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_2_4',0)",
                "VERTEX_NDOF_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_2_5',0)",
                "VERTEX_NDOF_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_2_6',0)",
                "VERTEX_NDOF_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_3_4',0)",
                "VERTEX_NDOF_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_3_5',0)",
                "VERTEX_NDOF_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_3_6',0)",
                "VERTEX_NDOF_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_4_5',0)",
                "VERTEX_NDOF_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_4_6',0)",
                "VERTEX_NDOF_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'VERTEX_NDOF_COMB_5_6',0)",

		"ETA_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_1_2',0)",
                "ETA_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_1_3',0)",
                "ETA_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_1_4',0)",
                "ETA_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_1_5',0)",
                "ETA_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_1_6',0)",
                "ETA_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_2_3',0)",
                "ETA_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_2_4',0)",
                "ETA_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_2_5',0)",
                "ETA_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_2_6',0)",
                "ETA_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_3_4',0)",
                "ETA_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_3_5',0)",
                "ETA_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_3_6',0)",
                "ETA_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_4_5',0)",
                "ETA_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_4_6',0)",
                "ETA_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'ETA_COMB_5_6',0)",


                "MCORR_OWNPV_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_1_2',0)",
                "MCORR_OWNPV_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_1_3',0)",
                "MCORR_OWNPV_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_1_4',0)",
                "MCORR_OWNPV_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_1_5',0)",
                "MCORR_OWNPV_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_1_6',0)",
                "MCORR_OWNPV_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_2_3',0)",
                "MCORR_OWNPV_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_2_4',0)",
                "MCORR_OWNPV_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_2_5',0)",
                "MCORR_OWNPV_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_2_6',0)",
                "MCORR_OWNPV_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_3_4',0)",
                "MCORR_OWNPV_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_3_5',0)",
                "MCORR_OWNPV_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_3_6',0)",
                "MCORR_OWNPV_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_4_5',0)",
                "MCORR_OWNPV_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_4_6',0)",
                "MCORR_OWNPV_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'MCORR_OWNPV_COMB_5_6',0)",


                "SUMPT_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_1_2',0)",
                "SUMPT_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_1_3',0)",
                "SUMPT_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_1_4',0)",
                "SUMPT_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_1_5',0)",
                "SUMPT_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_1_6',0)",
                "SUMPT_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_2_3',0)",
                "SUMPT_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_2_4',0)",
                "SUMPT_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_2_5',0)",
                "SUMPT_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_2_6',0)",
                "SUMPT_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_3_4',0)",
                "SUMPT_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_3_5',0)",
                "SUMPT_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_3_6',0)",
                "SUMPT_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_4_5',0)",
                "SUMPT_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_4_6',0)",
                "SUMPT_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'SUMPT_COMB_5_6',0)",

                "DIRA_OWNPV_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_1_2',0)",
                "DIRA_OWNPV_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_1_3',0)",
                "DIRA_OWNPV_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_1_4',0)",
                "DIRA_OWNPV_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_1_5',0)",
                "DIRA_OWNPV_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_1_6',0)",
                "DIRA_OWNPV_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_2_3',0)",
                "DIRA_OWNPV_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_2_4',0)",
                "DIRA_OWNPV_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_2_5',0)",
                "DIRA_OWNPV_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_2_6',0)",
                "DIRA_OWNPV_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_3_4',0)",
                "DIRA_OWNPV_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_3_5',0)",
                "DIRA_OWNPV_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_3_6',0)",
                "DIRA_OWNPV_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_4_5',0)",
                "DIRA_OWNPV_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_4_6',0)",
                "DIRA_OWNPV_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'DIRA_OWNPV_COMB_5_6',0)",


                "DOCA_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_1_2',0)",
                "DOCA_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_1_3',0)",
                "DOCA_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_1_4',0)",
                "DOCA_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_1_5',0)",
                "DOCA_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_1_6',0)",
                "DOCA_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_2_3',0)",
		"DOCA_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_2_4',0)",
                "DOCA_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_2_5',0)",
                "DOCA_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_2_6',0)",
                "DOCA_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_3_4',0)",
                "DOCA_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_3_5',0)",
                "DOCA_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_3_6',0)",
                "DOCA_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_4_5',0)",
                "DOCA_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_4_6',0)",
                "DOCA_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'DOCA_COMB_5_6',0)",

                "VDCHI2_OWNPV_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_1_2',0)",
                "VDCHI2_OWNPV_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_1_3',0)",
                "VDCHI2_OWNPV_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_1_4',0)",
                "VDCHI2_OWNPV_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_1_5',0)",
                "VDCHI2_OWNPV_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_1_6',0)",
                "VDCHI2_OWNPV_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_2_3',0)",
                "VDCHI2_OWNPV_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_2_4',0)",
                "VDCHI2_OWNPV_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_2_5',0)",
                "VDCHI2_OWNPV_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_2_6',0)",
                "VDCHI2_OWNPV_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_3_4',0)",
                "VDCHI2_OWNPV_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_3_5',0)",
                "VDCHI2_OWNPV_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_3_6',0)",
                "VDCHI2_OWNPV_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_4_5',0)",
                "VDCHI2_OWNPV_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_4_6',0)",
                "VDCHI2_OWNPV_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'VDCHI2_OWNPV_COMB_5_6',0)",


                "IPCHI2_OWNPV_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_1_2',0)",
                "IPCHI2_OWNPV_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_1_3',0)",
                "IPCHI2_OWNPV_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_1_4',0)",
                "IPCHI2_OWNPV_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_1_5',0)",
                "IPCHI2_OWNPV_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_1_6',0)",
                "IPCHI2_OWNPV_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_2_3',0)",
                "IPCHI2_OWNPV_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_2_4',0)",
                "IPCHI2_OWNPV_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_2_5',0)",
                "IPCHI2_OWNPV_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_2_6',0)",
                "IPCHI2_OWNPV_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_3_4',0)",
                "IPCHI2_OWNPV_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_3_5',0)",
                "IPCHI2_OWNPV_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_3_6',0)",
                "IPCHI2_OWNPV_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_4_5',0)",
                "IPCHI2_OWNPV_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_4_6',0)",
                "IPCHI2_OWNPV_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'IPCHI2_OWNPV_COMB_5_6',0)",

                "NLT_OWNPV_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_1_2',0)",
                "NLT_OWNPV_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_1_3',0)",
                "NLT_OWNPV_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_1_4',0)",
                "NLT_OWNPV_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_1_5',0)",
                "NLT_OWNPV_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_1_6',0)",
                "NLT_OWNPV_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_2_3',0)",
                "NLT_OWNPV_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_2_4',0)",
                "NLT_OWNPV_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_2_5',0)",
                "NLT_OWNPV_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_2_6',0)",
                "NLT_OWNPV_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_3_4',0)",
                "NLT_OWNPV_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_3_5',0)",
                "NLT_OWNPV_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_3_6',0)",
                "NLT_OWNPV_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_4_5',0)",
                "NLT_OWNPV_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_4_6',0)",
                "NLT_OWNPV_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'NLT_OWNPV_COMB_5_6',0)",

                "PT_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_1_2',0)",
                "PT_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_1_3',0)",
                "PT_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_1_4',0)",
                "PT_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_1_5',0)",
                "PT_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_1_6',0)",
                "PT_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_2_3',0)",
                "PT_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_2_4',0)",
                "PT_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_2_5',0)",
                "PT_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_2_6',0)",
                "PT_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_3_4',0)",
                "PT_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_3_5',0)",
                "PT_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_3_6',0)",
                "PT_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_4_5',0)",
                "PT_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_4_6',0)",
                "PT_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'PT_COMB_5_6',0)",

                "P_COMB_1_2" : "RELINFO('"+ relinfoutput+"', 'P_COMB_1_2',0)",
                "P_COMB_1_3" : "RELINFO('"+ relinfoutput+"', 'P_COMB_1_3',0)",
                "P_COMB_1_4" : "RELINFO('"+ relinfoutput+"', 'P_COMB_1_4',0)",
                "P_COMB_1_5" : "RELINFO('"+ relinfoutput+"', 'P_COMB_1_5',0)",
                "P_COMB_1_6" : "RELINFO('"+ relinfoutput+"', 'P_COMB_1_6',0)",
                "P_COMB_2_3" : "RELINFO('"+ relinfoutput+"', 'P_COMB_2_3',0)",
                "P_COMB_2_4" : "RELINFO('"+ relinfoutput+"', 'P_COMB_2_4',0)",
                "P_COMB_2_5" : "RELINFO('"+ relinfoutput+"', 'P_COMB_2_5',0)",
                "P_COMB_2_6" : "RELINFO('"+ relinfoutput+"', 'P_COMB_2_6',0)",
                "P_COMB_3_4" : "RELINFO('"+ relinfoutput+"', 'P_COMB_3_4',0)",
                "P_COMB_3_5" : "RELINFO('"+ relinfoutput+"', 'P_COMB_3_5',0)",
                "P_COMB_3_6" : "RELINFO('"+ relinfoutput+"', 'P_COMB_3_6',0)",
                "P_COMB_4_5" : "RELINFO('"+ relinfoutput+"', 'P_COMB_4_5',0)",
                "P_COMB_4_6" : "RELINFO('"+ relinfoutput+"', 'P_COMB_4_6',0)",
                "P_COMB_5_6" : "RELINFO('"+ relinfoutput+"', 'P_COMB_5_6',0)",

                "NCOMBINATIONS" : "RELINFO('"+ relinfoutput+"', 'NCOMBINATIONS', 0)"
                }



#This deletes stripping objects from previous stripping versions in the DST file                                     
from Configurables import EventNodeKiller
eventNodeKiller = EventNodeKiller('Stripkiller')
eventNodeKiller.Nodes = [ '/Event/AllStreams', '/Event/Strip' ]


