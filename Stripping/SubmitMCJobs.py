import sys
print (sys.argv)
import argparse

parser = argparse.ArgumentParser(description="Make my DaVinci job.")
parser.add_argument('--sample', choices=['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds','B_Lcpbarmunu','Lb_Lc2765munu','Lb_Lc2880munu'],
		                    help='which sample the script is run on', required=False)
parser.add_argument('--polarity', choices=['MagUp', 'MagDown'],
		                    help='Polarity of data-taking to run over', required=False)
parser.add_argument('--test', action='store_true',
				                    help='Run over one file locally')
parser.add_argument('--full', action='store_true',
				                    help='Run over full MC samples')

args = parser.parse_args()

sample = args.sample
polarity = args.polarity
test = args.test
fullsim = args.full
print(fullsim)

#print ('Running on : ', sample, polarity,'   performing test? ', test)

event = {'Lb_Lcmunu':'15874003',
	     'Lb_Lctaunu':'15874004',
	     'Lb_LcDs':'15894600',
	     'Lb_Lc2625taunu':'15576003',
	     'Lb_Lc2625munu':'15576002',
	     'Lb_Lc2593taunu':'15576005',
	     'Lb_Lc2593munu':'15576004',
	     'Lb_Lc2625Ds':'15896005',
	     'Lb_Lc2593Ds':'15896600',
         'B_Lcpbarmunu':'12575070',
         'Lb_Lc2765munu':'15876030',
         'Lb_Lc2880munu':'15876031'
	 }
samples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds']
samples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2593munu','Lb_Lc2593taunu','Lb_Lc2593Ds','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2625Ds','B_Lcpbarmunu','Lb_Lc2765munu','Lb_Lc2880munu']
#samples = ['Lb_Lc2880munu']
polarities = ['MagUp','MagDown']

def createjob(sample,polarity,bkPath,fullsim):
	print('Creating job: ',sample,' polarity: ',polarity)
	data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

	j = Job()
	j.name = sample+'_'+polarity
	j.comment = 'Wrong HLT1 MVA + true iso grandmaID information'
	myApp = GaudiExec()
	myApp.directory = "/home/hep/buonaura/DaVinciDev/DaVinciDev_v42r6p1/"
	j.application = myApp
	if (fullsim):
		j.application.options = ['RunStripping_OnMCfull.py']
	else:
		j.application.options = ['RunStripping_OnMC.py']
	j.application.platform = 'x86_64-centos7-gcc62-opt'
	j.backend = Dirac()
	j.outputfiles = [LocalFile('*.root')]
	j.inputfiles = [LocalFile('weights.xml')]
	if test==True:
		j.inputdata = data[0]
	else:
		j.inputdata = data
	#j.splitter = SplitByFiles(filesPerJob=2, ignoremissing = True)
	j.splitter = SplitByFiles(filesPerJob=1, ignoremissing = True)
	j.submit()
	return

if fullsim==False:
    for sample in samples:
        print(sample)
        for polarity in polarities:
            print(polarity)
            bkPath = '/MC/2016/Beam6500GeV-2016-'+polarity+'-TrackerOnly-Nu1.6-25ns-Pythia8/Sim09f/Reco16/Stripping28r1Filtered/'+event[sample]+'/LCTAUNU.SAFESTRIP.DST'
            if sample=='B_Lcpbarmunu' or sample=='Lb_Lc2765munu' or sample=='Lb_Lc2880munu':
                bkPath = '/MC/2016/Beam6500GeV-2016-'+polarity+'-TrackerOnly-Nu1.6-25ns-Pythia8/Sim09h/Reco16/Stripping28r1Filtered/'+event[sample]+'/LCTAUNU.SAFESTRIP.DST'
            print(bkPath)
            createjob(sample,polarity, bkPath,fullsim)
else:
	print('------- PROCESSING FULL SIMULATION FILES -----------')
	for sample in samples:
		print ('----- Submitting sample: ', sample)
		for polarity in polarities:
			print('---- POLARITY: ', polarity)
			bkPath = '/MC/2016/Beam6500GeV-2016-'+polarity+'-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28Filtered/'+event[sample]+'/LCTAUNU.SAFESTRIP.DST'
			createjob(sample,polarity, bkPath,fullsim)


