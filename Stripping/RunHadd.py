import sys, os
import argparse

parser = argparse.ArgumentParser(description="Add subjobs from DaVinci jobs.")
parser.add_argument('--datatype', choices=['MC','MCtrackeronly'])
args = parser.parse_args()

datatype = args.datatype

folder = {'MC':{'Lb_Lcmunu':{'MagUp':'200','MagDown':'201'},
                'Lb_Lctaunu':{'MagUp':'196','MagDown':'197'},
                'Lb_LcDs':{'MagUp':'206','MagDown':'207'},
                'Lb_Lc2625taunu':{'MagUp':'194','MagDown':'195'},
                'Lb_Lc2625munu':{'MagUp':'204','MagDown':'205'},
                'Lb_Lc2593taunu':{'MagUp':'198','MagDown':'199'},
                'Lb_Lc2593munu':{'MagUp':'202','MagDown':'203'}
                },
	   'MCtrackeronly':{'Lb_Lcmunu':{'MagUp':'277','MagDown':'278'},
	                   'Lb_Lctaunu':{'MagUp':'279','MagDown':'280'},
	                   'Lb_LcDs':{'MagUp':'295','MagDown':'296'},
	                   'Lb_Lc2625taunu':{'MagUp':'283','MagDown':'284'},
	                   'Lb_Lc2625munu':{'MagUp':'285','MagDown':'286'},
	                   'Lb_Lc2593taunu':{'MagUp':'287','MagDown':'288'},
	                   'Lb_Lc2593munu':{'MagUp':'289','MagDown':'290'},
	                   'Lb_Lc2625Ds':{'MagUp':'297','MagDown':'298'},
	                   'Lb_Lc2593Ds':{'MagUp':'299','MagDown':'300'}}
	}


datatypes = [datatype]
samples = {'MC':{'Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu'},
	   'MCtrackeronly':{'Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds'}}
polarities=['MagUp','MagDown']
fname={'MC':'tupleoutMC_trackeronly.root','Data':'tupleout.root','MCtrackeronly':'tupleoutMC_trackeronly.root'}

for dtype in datatypes:
    for sample in samples[dtype]:
	    #print(sample)
	    for polarity in polarities:
		    #print(polarity)
		    print('%s %s %s %s %s' %(folder[dtype][sample][polarity],sample, polarity, dtype,fname[dtype]))
           	    os.system('python haddTuples.py %s %s %s %s %s' %(folder[dtype][sample][polarity],sample, polarity, dtype,fname[dtype]))