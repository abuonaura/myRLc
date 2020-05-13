import sys, os
sys.path.append('L0Hadron_TOS/')
import emulate_L0HadronTOS
sys.path.append('HLT1/')
import emulate_HLT1_cuts_RDataframes
import emulateHLT1_twoTrack_RDataFrame
sys.path.append('L0Global_TIS/')
import emulate_L0GlobalTIS
from argparse import ArgumentParser

#inputFile = "/disk/lhcb_data2/ibezshyi/Lb_Lcmunu_MagUp.root"
inputFile = "/home/hep/ibezshyi/DaVinciDev_v45r1/DaVinciDev_v42r6p1/Lb_Lctaunu_MagDown.root"


parser = ArgumentParser()
parser.add_argument('-f',dest='inputFile', help="Input file if not default file", required=True, default=False)
parser.add_argument('--L0HTOS',dest='restartL0HTOS', help="Forces to reproduce L0Hadron TOS file", required=False, default=False, action='store_true')
parser.add_argument('--L0GTIS',dest='restartL0GTIS', help="Forces to reproduce L0Global TIS file", required=False, default=False, action='store_true')
parser.add_argument('--HLT1_1',dest='restartHLT1_1', help="Forces to reproduce HLT1 one track file", required=False, default=False, action='store_true')
parser.add_argument('--HLT1_2',dest='restartHLT1_2', help="Forces to reproduce HLT1 two track file", required=False, default=False, action='store_true')
options = parser.parse_args()

if options.inputFile:
    inputFile = options.inputFile
restartL0HTOS = options.restartL0HTOS
restartL0GTIS = options.restartL0GTIS
restartHLT1_1 = options.restartHLT1_1
restartHLT1_2 = options.restartHLT1_2
'''
print('Restarting L0 HadronTOS: ',restartL0HTOS)
print('Restarting L0 GlobalTIS: ',restartL0GTIS)
print('Restarting HLT1 One Track: ',restartHLT1_1)
print('Restarting HLT1 Two Track: ',restartHLT1_2)
'''

filedir = '/disk/lhcb_data2/RLcMuonic2016/MC/'
if restartL0HTOS==True:
    #print('Emulating LOHadronTOS')
    emulate_L0HadronTOS.main(inputFile)
if restartHLT1_1==True:
    #print('Emulating HLT1 One Track')
    emulate_HLT1_cuts_RDataframes.main(inputFile)
if restartHLT1_2==True:
    #print('Emulating HLT1 Two Track')
    emulateHLT1_twoTrack_RDataFrame.main(inputFile)
if restartL0GTIS==True:
    #print('Emulating L0GlobalTIS')
    emulate_L0GlobalTIS.main(inputFile)
    
