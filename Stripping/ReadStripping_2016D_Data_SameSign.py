#!/usr/bin/env gaudirun.py
# 
# Melody Ravonel Salzgeber
# Filter on various stripping line:
# -B2DMuNuX_Dp
# -B2DMuNuX_Ds
# -b2DsMuXB2DMuForTauMuLine
# 2017
# 
import GaudiKernel.SystemOfUnits as Units
from Gaudi.Configuration import *
from Configurables import DaVinci
from StrippingConf.StrippingStream import StrippingStream
from StrippingConf.Configuration import StrippingConf

MessageSvc().Format = "% F%60W%S%7W%R%T %0W%M"

from StrippingSettings.Stripping28.LineConfigDictionaries_Semileptonic import B2DMuForTauMu
from StrippingSelections.StrippingSL import StrippingB2DMuForTauMu
from GaudiConfUtils.ConfigurableGenerators import FilterDesktop, CombineParticles
from PhysSelPython.Wrappers import Selection
from StrippingConf.StrippingLine import StrippingLine
from StrippingUtils.Utils import LineBuilder
#track objs with specific mass hypothesis
from StandardParticles import StdLoosePions, StdLooseMuons, StdLooseKaons, StdLooseProtons, StdLooseElectrons
#StdNoPIDsPions, StdNoPIDsKaons: Pi and K assigned mass of Pi and K
from StandardParticles import StdNoPIDsPions, StdNoPIDsKaons,StdNoPIDsProtons,StdNoPIDsMuons
from GaudiKernel.SystemOfUnits import MeV, GeV, cm, mm



from Configurables import DecayTreeTuple, FitDecayTrees, TupleToolRecoStats, TupleToolTrigger, TupleToolSubMass
from Configurables import TupleToolTISTOS, CondDB, SelDSTWriter
from Configurables import TupleToolTrackInfo, TupleToolRICHPid, TupleToolGeometry, TupleToolPid
from Configurables import TupleToolANNPID
from Configurables import TupleToolSLTruth
from Configurables import TupleToolTauMuDiscrVarsLcMassConstraint
from DecayTreeTuple.Configuration import *

tupleB = DecayTreeTuple("tupleout")
tupleB.Inputs = ["/Event/Semileptonic/Phys/b2LcMuXB2DMuForTauMuLine/Particles"]
tupleB.Decay= "[B+ -> ^(Lambda_c+ -> ^p+ ^K- ^pi+) ^mu+]CC"
tupleB.ToolList =  [
      "TupleToolKinematic"
      , "TupleToolEventInfo"
      , "TupleToolRecoStats"
      #, "TupleToolGeometry"
      #,"TupleToolSLTruth"
] # Probably need to add many more Tools.

tupleB.addBranches ({
        "Lb": "[B+ -> (Lambda_c+ -> p+ K- pi+) mu+]CC",
        "Lc" : "[B+ -> ^Lambda_c+ mu+]CC",
         "p"  : "[B+ -> (Lambda_c+ -> ^p+ K- pi+) mu+]CC",
        "pi" : "[B+ -> (Lambda_c+ -> p+ K- ^pi+) mu+]CC",
        "K" : "[B+ -> (Lambda_c+ -> p+ ^K- pi+) mu+]CC",
        "mu" : "[B+ -> (Lambda_c+ -> p+ K- pi+) ^mu+]CC"
            })



l0_lines = [
  'L0HadronDecision'
  ,'L0MuonDecision'
  ,'L0ElectronDecision'
  ]

hlt1_lines = [
  'Hlt1TrackMVADecision'
  ,'Hlt1TwoTrackMVADecision'
  ]

hlt2_lines = [
  "Hlt2XcMuXForTauB2XcMuDecision",
  "Hlt2XcMuXForTauB2XcFakeMuDecision"
  ]


LoKi_All=tupleB.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_All")
LoKi_All.Variables = {
		'MINIPCHI2' : "MIPCHI2DV(PRIMARY)",
		'MINIP' : "MIPDV(PRIMARY)",
		'IPCHI2_OWNPV' : "BPVIPCHI2()",
		'IP_OWNPV' : "BPVIP()",
		'ghost' : "TRGHP",
		'ENDVERTEX_CHI2' : "VFASPF(VCHI2/VDOF)",
		'DIRA_OWNPV' : "BPVDIRA",
		'FDCHI2_OWNPV' : "BPVVDCHI2",
		'TRACK_CHI2' : "TRCHI2DOF"
		}

tt_taumudiscrvars_b = tupleB.Lb.addTupleTool("TupleToolTauMuDiscrVarsLcMassConstraint")

LoKi_B=tupleB.Lb.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_B")
LoKi_B.Variables = {
       'DOCAMAX' : "DOCAMAX",
       'TAU' : "BPVLTIME()",
       #'DIRA_OWNPV' : "BPVDIRA",
       #'FD_CHI2' : "BPVVDCHI2",
       #'ENDVERTEX_CHI2' : "VFASPF(VCHI2/VDOF)",
       'Corrected_Mass' : "BPVCORRM"
}

LoKi_B=tupleB.Lc.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_B")
LoKi_B.Variables = {
       'TAU' : "BPVLTIME()"
}

tupletistos = tupleB.Lb.addTupleTool("TupleToolTISTOS")
tupletistos.VerboseL0 = True
tupletistos.VerboseHlt1 = False
tupletistos.VerboseHlt2 = True
tupletistos.TriggerList = l0_lines+hlt2_lines

tupletistos1= tupleB.Lc.addTupleTool("TupleToolTISTOS")
tupletistos1.VerboseL0 = True
tupletistos1.VerboseHlt1 = True
tupletistos1.VerboseHlt2 = False
tupletistos1.TriggerList = l0_lines + hlt1_lines 

LoKi_muplus=tupleB.mu.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_muplus")
LoKi_muplus.Variables = {
       'PIDmu' : "PIDmu",
       'NNmu' : "PPINFO(PROBNNmu)"
}

#do similar for kaon,proton,pion for ProbNNpi/K/p

#****
#To configure the TupleTool:
#****
from Configurables import TupleToolApplyIsolation
tupleB.Lb.addTool(TupleToolApplyIsolation, name="TupleToolApplyIsolation")
tupleB.Lb.TupleToolApplyIsolation.WeightsFile="weights.xml"
tupleB.Lb.ToolList+=["TupleToolApplyIsolation/TupleToolApplyIsolation"]
#****

#TupleTool for Lc lifetime
tupleB.Lc.ToolList += [ "TupleToolPropertime" ]


#Add TupleTool with DLL                                                                                                                                                                                                                     
tupleB.mu.addTool(TupleToolPid,name='TupleToolPid')
tupleB.mu.ToolList+=['TupleToolPid']

tupleB.pi.addTool(TupleToolPid,name='TupleToolPid')
tupleB.pi.ToolList+=['TupleToolPid']

tupleB.K.addTool(TupleToolPid,name='TupleToolPid')
tupleB.K.ToolList+=['TupleToolPid']

tupleB.p.addTool(TupleToolPid,name='TupleToolPid')
tupleB.p.ToolList+=['TupleToolPid']


#Add TupleTool with Yandex PID
'''
tupleB.mu.addTool(TupleToolANNPID,name='TupleToolANNPID')
tupleB.mu.ToolList+=['TupleToolANNPID']
tupleB.mu.TupleToolANNPID.ANNPIDTunes = ['MC15TuneV1','MC15TuneDNNV1','MC15TuneFLAT4dV1','MC15TuneCatBoostV1']

tupleB.pi.addTool(TupleToolANNPID,name='TupleToolANNPID')
tupleB.pi.ToolList+=['TupleToolANNPID']
tupleB.pi.TupleToolANNPID.ANNPIDTunes = ['MC15TuneV1','MC15TuneDNNV1','MC15TuneFLAT4dV1','MC15TuneCatBoostV1']

tupleB.K.addTool(TupleToolANNPID,name='TupleToolANNPID')
tupleB.K.ToolList+=['TupleToolANNPID']
tupleB.K.TupleToolANNPID.ANNPIDTunes = ['MC15TuneV1','MC15TuneDNNV1','MC15TuneFLAT4dV1','MC15TuneCatBoostV1']

tupleB.p.addTool(TupleToolANNPID,name='TupleToolANNPID')
tupleB.p.ToolList+=['TupleToolANNPID']
tupleB.p.TupleToolANNPID.ANNPIDTunes = ['MC15TuneV1','MC15TuneDNNV1','MC15TuneFLAT4dV1','MC15TuneCatBoostV1']
'''

#To produce the Velo particles input container (the gain in performance from adding velo tracks is actually pretty negligible...)
#*********************************

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





#This deletes stripping objects from previous stripping versions in the DST file
#from Configurables import EventNodeKiller
#eventNodeKiller = EventNodeKiller('Stripkiller')
#eventNodeKiller.Nodes = [ '/Event/AllStreams', '/Event/Strip' ]

DaVinci().appendToMainSequence( [ tupleB] )
DaVinci().HistogramFile = 'DV_stripping_histos.root'
DaVinci().EvtMax = -1
#DaVinci().EvtMax = 2000
DaVinci().PrintFreq = 1000

#DaVinci().ProductionType = "Stripping"
DaVinci().DataType  = "2016"
DaVinci().InputType = "DST"

# database
#what are database condition of the file you're stripping on
#DaVinci().DDDBtag   = "dddb-20150724"
#DaVinci().CondDBtag = 'sim-20161124-2-vc-md100'

DaVinci().Simulation = False
DaVinci().Lumi =  not DaVinci().Simulation

DaVinci().TupleFile = "tupleout.root"

# change the column size of timing table
from Configurables import TimingAuditor, SequencerTimerTool
TimingAuditor().addTool(SequencerTimerTool,name="TIMER")
TimingAuditor().TIMER.NameSize = 60


#DaVinci().Input=['PFN:/home/hep/buonaura/Analysis/RLc/Datasets/Data/00069529_00003917_1.semileptonic.dst']  



#from GaudiConf import IOHelper

#IOHelper().inputFiles([
#    '/eos/lhcb/user/i/ibezshyi/ganga/00059560_00000019_1.semileptonic.dst'
#], clear=True)



