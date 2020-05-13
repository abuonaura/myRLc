import ROOT as r
import collections
import sys
from ReduceTree import *
from CreateFinalWorkingSamples import *

outf = open('Efficiencies.txt','w')

datatype = ['MC','Data','FakeMu','FakeMuSS','DataSS']
polarities = ['MagUp','MagDown']
mcsamples = ['LcTau','LcMu','LcDs','Lc2593Mu','Lc2593Tau','Lc2625Mu','Lc2625Tau']

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
startfile = {'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp.root','MagDown':datadir+'Data/Lb_Data_MagDown.root'},
        'DataSS':{'MagUp':datadir+'Data/Lb_DataSS_MagUp.root','MagDown':datadir+'Data/Lb_DataSS_MagDown.root'},
        'FakeMu':{'MagUp':datadir+'Data/Lb_FakeMu_MagUp.root',
                  'MagDown':datadir+'Data/Lb_FakeMu_MagDown.root'},
        'FakeMuSS':{'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp.root',
                  'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown.root'},
        'MC':{'LcTau':datadir+'MC/Lb_Lctaunu_PID.root',
              'LcMu':datadir+'MC/Lb_Lcmunu_PID.root',
              'LcDs':datadir+'MC/Lb_LcDs_PID.root',
              'Lc2593Mu':datadir+'MC/Lb_Lc2593munu_PID.root',
              'Lc2593Tau':datadir+'MC/Lb_Lc2593taunu_PID.root',
              'Lc2625Mu':datadir+'MC/Lb_Lc2625munu_PID.root',
              'Lc2625Tau':datadir+'MC/Lb_Lc2625taunu_PID.root'}}
reducedfile = {'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced.root','MagDown':datadir+'Data/Lb_Data_MagDown_reduced.root'},
        'DataSS':{'MagUp':datadir+'Data/Lb_DataSS_MagUp_reduced.root','MagDown':datadir+'Data/Lb_DataSS_MagDown_reduced.root'},
        'FakeMu':{'MagUp':datadir+'Data/Lb_FakeMu_MagUp_reduced.root',
                  'MagDown':datadir+'Data/Lb_FakeMu_MagDown_reduced.root'},
        'FakeMuSS':{'MagUp':datadir+'Data/Lb_FakeMuSS_MagUp_reduced.root',
                  'MagDown':datadir+'Data/Lb_FakeMuSS_MagDown_reduced.root'},
        'MC':{'LcTau':datadir+'MC/Lb_Lctaunu_PID_reduced.root',
              'LcMu':datadir+'MC/Lb_Lcmunu_PID_reduced.root',
              'LcDs':datadir+'MC/Lb_LcDs_PID_reduced.root',
              'Lc2593Mu':datadir+'MC/Lb_Lc2593munu_PID_reduced.root',
              'Lc2593Tau':datadir+'MC/Lb_Lc2593taunu_PID_reduced.root',
              'Lc2625Mu':datadir+'MC/Lb_Lc2625munu_PID_reduced.root',
              'Lc2625Tau':datadir+'MC/Lb_Lc2625taunu_PID_reduced.root',
              }}


totEvents=collections.defaultdict(dict)
triggerEvents=collections.defaultdict(dict)
LcMcutEvents=collections.defaultdict(dict)
MCpidEvents=collections.defaultdict(dict)
BDTEvents=collections.defaultdict(dict)
pProbNNEvents=collections.defaultdict(dict)
mucutEvents=collections.defaultdict(dict)

triggerEff = collections.defaultdict(dict)
LcMEff = collections.defaultdict(dict)
MCpidEff = collections.defaultdict(dict)
BDTEff = collections.defaultdict(dict)
pProbNNEff=collections.defaultdict(dict)
mucutEff=collections.defaultdict(dict)

cMass = r.TCut("Lc_M>2230 && Lc_M<2330")
BDTcut = GetBDTcut(0.65)
pProbNNcut = GetpProbNNcut()
mucut = GetmuPIDcut()


for dtype in datatype:
    triggercut = GetTriggerCut(dtype)
    pidcut = GetStrippingPIDcut()
    print('',file = outf)
    print('>>>>>>>>> '+dtype+' files <<<<<<<<<<',file=outf)
    if dtype!='MC':
        for polarity in polarities:
            print('',file=outf)
            print('   ---------  Polarity: '+ polarity+' ----------',file=outf)
            print('',file=outf)
            f = r.TFile(startfile[dtype][polarity],'READ')
            t = f.Get('tupleout/DecayTree')
            totEvents[dtype][polarity]= t.GetEntries()
            print(' Total Evts:         ',totEvents[dtype][polarity],file=outf)
            triggerEvents[dtype][polarity] = int(t.Draw('Lb_M',triggercut,'goff'))
            triggerEff[dtype][polarity] = triggerEvents[dtype][polarity]*1./totEvents[dtype][polarity]
            LcMcutEvents[dtype][polarity] = int(t.Draw('Lb_M',triggercut.GetTitle()+'&&'+cMass.GetTitle(),'goff'))
            LcMEff[dtype][polarity] = LcMcutEvents[dtype][polarity]*1./triggerEvents[dtype][polarity]
            print(' Trigger Cuts:       ',triggerEvents[dtype][polarity], '      eff: {:.4f}'.format(triggerEff[dtype][polarity]),file=outf)
            print(' Lc Mass Cut:        ',LcMcutEvents[dtype][polarity], '      eff: {:.4f}'.format(LcMEff[dtype][polarity]),file=outf)

            f.Close()

            f = r.TFile(reducedfile[dtype][polarity],'READ')
            t = f.Get('DecayTree')
            BDTEvents[dtype][polarity] = int(t.Draw('Lb_M',BDTcut,'goff'))
            pProbNNEvents[dtype][polarity] = int(t.Draw('Lb_M',BDTcut.GetTitle()+'&&'+pProbNNcut.GetTitle(),'goff'))
            BDTEff[dtype][polarity] = BDTEvents[dtype][polarity]*1./LcMcutEvents[dtype][polarity]
            pProbNNEff[dtype][polarity] = pProbNNEvents[dtype][polarity]*1./BDTEvents[dtype][polarity]
            print(' BDT Cut:            ',BDTEvents[dtype][polarity], '      eff: {:.4f}'.format(BDTEff[dtype][polarity]),file=outf)
            print(' pProbNN Cut:        ',pProbNNEvents[dtype][polarity], '      eff: {:.4f}'.format(pProbNNEff[dtype][polarity]),file=outf)

            if dtype=='Data' or dtype =='DataSS':
                mucutEvents[dtype][polarity] = int(t.Draw('Lb_M',BDTcut.GetTitle()+'&&'+pProbNNcut.GetTitle()+'&&'+mucut.GetTitle(),'goff'))
                mucutEff[dtype][polarity] = mucutEvents[dtype][polarity]*1./pProbNNEvents[dtype][polarity]
                print(' muPID Cut:          ',mucutEvents[dtype][polarity], '      eff: {:.4f}' .format(mucutEff[dtype][polarity]),file=outf)
                
    if dtype=='MC':
        for sample in mcsamples:
            print('',file=outf)
            print('   ---------  Sample: '+ sample +' ----------',file=outf)
            print('',file=outf)
            
            f = r.TFile(startfile[dtype][sample],'READ')
            t = f.Get('DecayTree')
            totEvents[dtype][sample]= t.GetEntries()
            triggerEvents[dtype][sample] = int(t.Draw('Lb_M',triggercut,'goff'))
            LcMcutEvents[dtype][sample] = int(t.Draw('Lb_M',triggercut.GetTitle()+'&&'+cMass.GetTitle(),'goff'))
            MCpidEvents[dtype][sample] = int(t.Draw('Lb_M',triggercut.GetTitle()+'&&'+cMass.GetTitle()+'&&'+pidcut.GetTitle(),'goff'))
            
            triggerEff[dtype][sample] = triggerEvents[dtype][sample]*1./totEvents[dtype][sample]
            LcMEff[dtype][sample] = LcMcutEvents[dtype][sample]*1./triggerEvents[dtype][sample]
            MCpidEff[dtype][sample] = MCpidEvents[dtype][sample]*1./LcMcutEvents[dtype][sample]
            
            f.Close()

            f = r.TFile(reducedfile[dtype][sample],'READ')
            t = f.Get('DecayTree')
            BDTEvents[dtype][sample] = int(t.Draw('Lb_M',BDTcut,'goff'))
            BDTEff[dtype][sample] = BDTEvents[dtype][sample]*1./LcMcutEvents[dtype][sample]
            print(' Total Evts:         ',totEvents[dtype][sample],file=outf)
            print(' Trigger Cuts:       ',triggerEvents[dtype][sample], '      eff: {:.4f}'.format(triggerEff[dtype][sample]),file=outf)
            print(' Lc Mass Cut:        ',LcMcutEvents[dtype][sample], '      eff: {:.4f}'.format(LcMEff[dtype][sample]),file=outf)
            print(' Stripping PID Cut:  ', MCpidEvents[dtype][sample], '      eff: {:.4f}'.format(MCpidEff[dtype][sample]),file=outf)
            print(' BDT Cut:            ',BDTEvents[dtype][sample], '      eff: {:.4f}'.format(BDTEff[dtype][sample]),file=outf)

outf.close()
