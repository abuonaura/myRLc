import os
import sys
import ROOT as r
import numpy as np
import os, sys, glob


#!!!! @nn !!!!#
#rescale nTracks in input ntuple

#samples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu']
#samples = ['Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu','Lb_Lc2593taunu','Lb_Lc2593munu']
samples = ['Lb_Lc2625Ds']
datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MC/'
#polarities = ['MagUp','MagDown']
polarities = ['MagDown']

for sample in samples:
    for polarity in polarities:
	startfile = datadir+sample+'_'+polarity+'.root'
	print('>>>>   Adding PID to %s with polarity %s' %(sample,polarity))
	print(startfile)
	
	modfile = startfile[:-5]+'_rescaled.root'
	pidfile = startfile[:-5]+'_PID.root'
	print(modfile, pidfile)


	print('... Rescaling ntracks ...')

	f = r.TFile(startfile,'READ')
	t = f.Get('tupleout/DecayTree')


	of = r.TFile(modfile,'RECREATE')
	ot = t.CloneTree(0)

	nTracks = np.zeros(1, dtype=int)
	ot.SetBranchAddress("nTracks",nTracks)

	for i in range(t.GetEntries()):
	    t.GetEntry(i)
	    nTracks[0] = t.nTracks*1.16
	    ot.Fill()


	ot.Write()
	of.Close()
	f.Close()

	print('... End of rescaling... ')
	print('... Adding PIDs ...' )

	## START OF CONFIG
	# Read comments and check vars
	# at least until end of config section

	# List of input ROOT files with MC ntuples. Format:
	#   (inputfile, outputfile, dataset)
	files = [
	  (modfile, pidfile, "MagDown_2016"),
	]

	# Name of the input tree
	# Could also include ROOT directory, e.g. "Dir/Ntuple"
	input_tree = "DecayTree"

	# Postfixes of the Pt, Eta and Ntracks variables (ntuple variable name w/o branch name)
	# e.g. if the ntuple contains "pion_PT", it should be just "PT"
	ptvar  = "PT"
	#etavar = "eta"
	#pvar   = None
	## Could use P variable instead of eta
	etavar = None
	pvar   = "P"

	ntrvar = "nTracks"  # This should correspond to the number of "Best tracks", not "Long tracks"!

	seed = None   # No initial seed
	# seed = 1    # Alternatively, could set initial random seed

	# Dictionary of tracks with their PID variables, in the form {branch name}:{pidvars}
	# For each track branch name, {pidvars} is a dictionary in the form {ntuple variable}:{pid config},
	#   where
	#     {ntuple variable} is the name of the corresponding ntuple PID variable without branch name,
	#   and
	#     {pid_config} is the string describing the PID configuration.
	# Run PIDCorr.py without arguments to get the full list of PID configs

	tracks = {
	  'K'   : {
		    "PIDK"  : "K_CombDLLK_Brunel",
		    "PIDmu"  : "K_CombDLLmu_Brunel",
		    "PIDp"  : "K_CombDLLp_Brunel",
		    "ProbNNk" : "K_MC15TuneV1_ProbNNK_Brunel",
		    "ProbNNmu" : "K_MC15TuneV1_ProbNNmu_Brunel",
		    "ProbNNp" : "K_MC15TuneV1_ProbNNp_Brunel",
		    "ProbNNpi" : "K_MC15TuneV1_ProbNNpi_Brunel",
		   },
	  'mu' : {
		    "PIDK"  : "mu_CombDLLK_Brunel",
		    "PIDmu"  : "mu_CombDLLmu_Brunel",
		    "ProbNNk" : "mu_MC15TuneV1_ProbNNK_Brunel",
		    "ProbNNpi" : "mu_MC15TuneV1_ProbNNpi_Brunel",
		    "ProbNNmu" : "mu_MC15TuneV1_ProbNNmu_Brunel",
		   },
	  'p'  : {
		    "PIDK"  : "p_CombDLLK_Brunel",
		    "PIDp"  : "p_CombDLLp_Brunel",
		    "ProbNNk" : "p_MC15TuneV1_ProbNNK_Brunel",
		    "ProbNNpi" : "p_MC15TuneV1_ProbNNpi_Brunel",
		    "ProbNNp" : "p_MC15TuneV1_ProbNNp_Brunel",
		   },
	  'pi'  : {
		    "PIDK" : "pi_CombDLLK_Brunel",
		    "PIDmu" : "pi_CombDLLmu_Brunel",
		    "PIDp" : "pi_CombDLLp_Brunel",
		    "ProbNNk" : "pi_MC15TuneV1_ProbNNK_Brunel",
		    "ProbNNmu" : " pi_MC15TuneV1_ProbNNmu_Brunel",
		    "ProbNNp" : "pi_MC15TuneV1_ProbNNp_Brunel",
		    "ProbNNpi" : "pi_MC15TuneV1_ProbNNpi_Brunel",
		   },
	}


	# IF ON LXPLUS: if /tmp exists and is accessible, use for faster processing
	# IF NOT: use /tmp if you have enough RAM
	# temp_folder = '/tmp'
	# ELSE: use current folder
	temp_folder = '.'

	## END OF CONFIG


	# make sure we don't overwrite local files and prefix them with random strings
	import string
	import random
	rand_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))  # get 10 random chars for temp_file prefix

	temp_file_prefix = temp_folder + '/' + rand_string  # prefix temp files with folder and unique ID

	output_tree = input_tree.split("/")[-1]
	treename = input_tree

	for input_file, output_file, dataset in files :
	  tmpinfile = input_file
	  tmpoutfile = "%s_tmp1.root" % temp_file_prefix
	  for track, subst in tracks.iteritems() :
	    for var, config in subst.iteritems() :
	      command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py"
	      command += " -m %s_%s" % (track, ptvar)
	      if etavar:
                  command += " -e %s_%s" % (track, etavar)
	      elif pvar:
                  command += " -q %s_%s" % (track, pvar)
	      else:
                  print('Specify either ETA or P branch name per track')
                  sys.exit(1)
	      command += " -n %s" % ntrvar
	      command += " -t %s" % treename
	      command += " -p %s_%s_corr" % (track, var)
	      command += " -c %s" % config
	      command += " -d %s" % dataset
	      command += " -i %s" % tmpinfile
	      command += " -o %s" % tmpoutfile
	      if seed :
                  command += " -s %d" % seed

	      treename = output_tree
	      tmpinfile = tmpoutfile
	      if 'tmp1' in tmpoutfile:
                  tmpoutfile = tmpoutfile.replace('tmp1', 'tmp2')
	      else :
                  tmpoutfile = tmpoutfile.replace('tmp2', 'tmp1')

	      print(command)
	      os.system(command)

	  if "root://" in output_file:
	    print("xrdcp %s %s" % (tmpinfile, output_file))
	    os.system("xrdcp %s %s" % (tmpinfile, output_file))
	  else:
	    print("mv %s %s" % (tmpinfile, output_file))
	    os.system("mv %s %s" % (tmpinfile, output_file))
