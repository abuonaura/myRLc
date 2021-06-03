import sys, os, glob
import argparse
import ROOT as r

parser = argparse.ArgumentParser(description="Add files produced by Ganga")
#Arguments without the '--' are mandatory
parser.add_argument('jobID',help = 'folder number')
parser.add_argument('jobname', choices=['Lb_Lcmunu', 'Lb_Lctaunu','Lb_LcDs','Lb_Lc2625taunu','Lb_Lc2625munu', 'Lb_Lc2593taunu', 'Lb_Lc2593munu','Lb_Lc2625Ds','Lb_Lc2593Ds','B_Lcpbarmunu','Lb_Lc2765munu','Lb_Lc2880munu','Lb_Data','Lb_FakeMu','Lb_DataSS','Lb_FakeMuSS'], help = 'name of the sample')
parser.add_argument('polarity',choices=['MagUp','MagDown'], help = 'sample polarity')
parser.add_argument('dtype',choices=['Data','MC','MCtrackeronly'], help = 'type: data or MC files')
parser.add_argument('filename', help = 'name of the tree')

args = parser.parse_args()

jobID = args.jobID
jobname = args.jobname
polarity = args.polarity
dtype = args.dtype
fname = args.filename


print('%s %s %s %s %s'%(jobID, jobname, polarity, dtype, fname))

#gangadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/'
gangadir = '/disk/lhcb_data2/buonaura/gangadir/workspace/buonaura/LocalXML/'

files2Hadd=['']


#Retrieve the list of files to hadd
infiles = glob.glob(gangadir + '%s/*/output/'%jobID+fname)
                        
nsubjobs = len(infiles)
print(infiles, nsubjobs)
files2Hadd=''

infiles2add=[]

#Verify that the input files are not zombie
for ifile in infiles:
    f= r.TFile(ifile, 'READ')
    if f:
        if not f.IsZombie():
            infiles2add.append(ifile)
    f.Close()
nfiles2add = len(infiles2add)
files2Hadd += ' '.join(infiles2add)
print ('\nJob %s: Found %s outputs in %s subjobs.' % (jobID, nfiles2add, nsubjobs))

print('.... Hadding files ...')

#Create name output file
if dtype=='MC':
	#outfile = '/disk/lhcb_data2/RLcMuonic2016/'+dtype+'_full_new/'+jobname+'_'+polarity+'.root'
	outfile = '/disk/lhcb_data2/RLcMuonic2016/'+dtype+'_full_trueTrigger/'+jobname+'_'+polarity+'.root'
if dtype=='MCtrackeronly':
	outfile = '/disk/lhcb_data2/RLcMuonic2016/MC_TrackerOnly/'+jobname+'_'+polarity+'.root'
print(outfile)

logfile='HaddLog_'+jobname+'_'+polarity+'.txt'
print(' >>>> Creating %s file with polarity %s' %(jobname,polarity))
if os.path.isfile(outfile):
    print ('File already there... removing old file')
    os.system('rm %s' %outfile)
os.system('hadd -f6 %s %s >> %s' %(outfile, files2Hadd,logfile))


