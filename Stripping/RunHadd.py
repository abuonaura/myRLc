import sys, os
import argparse

parser = argparse.ArgumentParser(description="Add subjobs from DaVinci jobs.")
parser.add_argument('--datatype', choices=['MC','MCtrackeronly'])
args = parser.parse_args()

datatype = args.datatype

folder = {'MC':{'Lb_Lcmunu':{'MagUp':'361','MagDown':'362'},
                'Lb_Lctaunu':{'MagUp':'363','MagDown':'364'},
                'Lb_LcDs':{'MagUp':'365','MagDown':'366'},
                'Lb_Lc2593munu':{'MagUp':'367','MagDown':'368'},
                'Lb_Lc2593taunu':{'MagUp':'369','MagDown':'370'},
                'Lb_Lc2593Ds':{'MagUp':'371','MagDown':'372'},
                'Lb_Lc2625munu':{'MagUp':'373','MagDown':'374'},
                'Lb_Lc2625taunu':{'MagUp':'375','MagDown':'376'},
                'Lb_Lc2625Ds':{'MagUp':'377','MagDown':'378'}
                },
	   'MCtrackeronly':{'Lb_Lcmunu':{'MagUp':'343','MagDown':'344'},
	                   'Lb_Lctaunu':{'MagUp':'345','MagDown':'346'},
	                   'Lb_LcDs':{'MagUp':'347','MagDown':'348'},
	                   'Lb_Lc2593munu':{'MagUp':'349','MagDown':'350'},
	                   'Lb_Lc2593taunu':{'MagUp':'351','MagDown':'352'},
	                   'Lb_Lc2593Ds':{'MagUp':'353','MagDown':'354'},
	                   'Lb_Lc2625munu':{'MagUp':'355','MagDown':'356'},
	                   'Lb_Lc2625taunu':{'MagUp':'357','MagDown':'358'},
	                   'Lb_Lc2625Ds':{'MagUp':'359','MagDown':'360'}}
	}


datatypes = [datatype]
samples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds']
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
