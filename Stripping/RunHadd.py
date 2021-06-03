import sys, os
import argparse

parser = argparse.ArgumentParser(description="Add subjobs from DaVinci jobs.")
parser.add_argument('--datatype', choices=['MC','MCtrackeronly'])
args = parser.parse_args()

datatype = args.datatype

folder = {'MC':{'Lb_Lcmunu':{'MagUp':'435','MagDown':'436'},
                'Lb_Lctaunu':{'MagUp':'437','MagDown':'438'},
                'Lb_LcDs':{'MagUp':'439','MagDown':'440'},
                'Lb_Lc2593munu':{'MagUp':'441','MagDown':'442'},
                'Lb_Lc2593taunu':{'MagUp':'443','MagDown':'444'},
                'Lb_Lc2593Ds':{'MagUp':'445','MagDown':'446'},
                'Lb_Lc2625munu':{'MagUp':'447','MagDown':'448'},
                'Lb_Lc2625taunu':{'MagUp':'449','MagDown':'450'},
                'Lb_Lc2625Ds':{'MagUp':'451','MagDown':'452'}
                },
	   'MCtrackeronly':{'Lb_Lcmunu':    {'MagUp':'0','MagDown':'1'},
	                   'Lb_Lctaunu':    {'MagUp':'2','MagDown':'3'},
	                   'Lb_LcDs':       {'MagUp':'4','MagDown':'5'},
	                   'Lb_Lc2593munu': {'MagUp':'6','MagDown':'7'},
	                   'Lb_Lc2593taunu':{'MagUp':'8','MagDown':'9'},
	                   'Lb_Lc2593Ds':   {'MagUp':'10','MagDown':'11'},
	                   'Lb_Lc2625munu': {'MagUp':'12','MagDown':'13'},
	                   'Lb_Lc2625taunu':{'MagUp':'14','MagDown':'15'},
	                   'Lb_Lc2625Ds':   {'MagUp':'16','MagDown':'17'},
					   'B_Lcpbarmunu':  {'MagUp':'18','MagDown':'19'},
					   'Lb_Lc2765munu': {'MagUp':'20','MagDown':'21'},
					   'Lb_Lc2880munu': {'MagUp':'22','MagDown':'23'}
					   }
	}


datatypes = [datatype]
samples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds']
if datatype=='MCtrackeronly':
	samples += ['B_Lcpbarmunu','Lb_Lc2765munu','Lb_Lc2880munu']
polarities=['MagUp','MagDown']
fname={'MC':'tupleoutMC.root','Data':'tupleout.root','MCtrackeronly':'tupleoutMC.root'}
#fname = {'MC':'LbLcDs_fullsim.root'}

for dtype in datatypes:
	for sample in samples:
		print(sample)
		for polarity in polarities:
			print(polarity)
			print('%s %s %s %s %s' %(folder[dtype][sample][polarity],sample, polarity, dtype,fname[dtype]))
			os.system('python haddTuples.py %s %s %s %s %s' %(folder[dtype][sample][polarity],sample, polarity, dtype,fname[dtype]))
