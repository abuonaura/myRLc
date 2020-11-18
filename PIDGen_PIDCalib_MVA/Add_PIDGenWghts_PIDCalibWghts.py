import os
from ROOT import TFile, TTree, TH3F, TH1, kFALSE
import shutil
from root_pandas import read_root

magtypes = ['MagUp', 'MagDown']
years    = ['2011' , '2012', '2016', '2017', '2018']
trcks    = ['pi', 'K', 'p', 'mu']
leafs    = ['PT', 'P', 'nTracks']

perfhistname = {}
perfhistname['K']  = "K_DLLK>4_All"
perfhistname['Pi'] = "Pi_DLLK<2_All"
perfhistname['P']  = "P_DLLp>0_All"
perfhistname['Mu'] = "Mu_DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1&&MC15TuneV1_ProbNNghost<0.2_All"
trcks_PIDCalib     = ['Pi', 'K', 'P', 'Mu']
leafs_PIDCalib     = ['Brunel_P', 'Brunel_PT', 'nTracks_Brunel']

def Exceptions(infname, year, magtype, tmpdir, UraniaDir):
    if 'root://eoslhcb.cern.ch/' in infname:
        print('Using file in EOS space, check manually the file exists?')
    else:
        if not os.path.isfile(infname):
            raise Exception(infname+' does not exist!')

    if year not in years: 
        raise Exception(year+'not in '+years)

    if magtype not in magtypes: 
        raise Exception(magtype+'not in '+magtypes)

    if not os.path.isdir(tmpdir):
        raise Exception(tmpdir+' directory does not exist!')

    if not os.path.isdir(UraniaDir):
        raise Exception(UraniaDir+' directory does not exist!')

def make_tmpdir(tmpdir):
    if not tmpdir.endswith('/'): tmpdir += '/'
    tmpdir += 'tmp/'
    if not os.path.isdir(tmpdir): os.mkdir(tmpdir)
    return tmpdir

def make_tmpfile(infname, intreename, tmpinfname, varstokeep, nentries_to_read = 1000000000):
    print('Taking as input Tfile', infname,'with ttreename', intreename,'and storing tmp file', varstokeep)
    print('Only these vars are kept in the tmpfile', tmpinfname)
    ftfile = TFile(infname,'READ')
    ttree  = ftfile.Get(intreename)
    ttree.SetBranchStatus('*',0)
    for vartokeep in varstokeep: ttree.SetBranchStatus(vartokeep,1)
    tmpftfile = TFile(tmpinfname,'recreate')
    tmpttree  = ttree.CopyTree('', '', nentries_to_read)
    tmpttree.Write()
    tmpftfile.Close()
    ftfile.Close()

def AddPIDGenWeights(infname, intreename, outfname, magtype, UraniaDir, year = '2016', tmpdir = '/tmp/', nentries_to_read = 1000000000):
    """
    Function:
    ---------
    Adds PIDGen variables
    Proton: p_PROBNNp_corr
    Kaon  : K_PROBNNK_corr
    Pion  : pi_PROBNNpi_corr

    Arguments:
    ---------
    infname:    full path to input file tfile including it's name
    intreename: ttree name in the input file
    outfname:   full path to the output file including it's name where it will be written
    magtype :   MagDown or MagUp
    year :      2016
    tmpdir :    Currently default is set to /tmp/.The program makes a tmp directory inside tmpdir i.e. /tmp/tmp. It also later deltes it.
    UraniaDir:  Path to the Urania directory that holds the run bash script
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000

    Output:
    ------
    Stores a root file with the same TTree name as the input. It only contains eventNumber, runNumber and PIDGen vars.
    """
    Exceptions(infname, year, magtype, tmpdir, UraniaDir)

    print('Making tmp directory to temporarily store files in', tmpdir)
    tmpdir = make_tmpdir(tmpdir)
    if not UraniaDir.endswith('/'): UraniaDir += '/'

    print('Using the options', infname, intreename, outfname, magtype, year, tmpdir, UraniaDir)

    tmpfinfname = tmpdir+(infname.split('/')[-1]).replace('.root', '_pidgentmp.root')
    leafskeep   = [trck+'_'+lf for trck in trcks for lf in leafs[:2]] + ['runNumber', 'eventNumber', leafs[2]]
    make_tmpfile(infname, intreename, tmpfinfname, leafskeep, nentries_to_read)

    files       = [(tmpfinfname, outfname, magtype+"_"+year)]
    output_tree = intreename.split("/")[-1]

    ptvar      = leafs[0]
    pvar       = leafs[1]
    ntrvar     = leafs[2]
    etavar     = None

    seed = None   #seed = 1 #Alternatively, could set initial random seed

    tracks = {}
    tracks['pi']= {"PROBNNpi"   : "pi_PIDKlt2_MC15TuneV1_ProbNNpi_Brunel"}
    tracks['K'] = {"PROBNNK"    : "K_PIDKgt4_MC15TuneV1_ProbNNK_Brunel"}
    tracks['p'] = {"PROBNNp"    : "p_LbLcPi_PIDpgt0_MC15TuneV1_ProbNNp_Brunel"}

    for input_file, output_file, dataset in files : 
        #treename   = intreename
        treename   = "DecayTree"
        tmpinfile  = input_file
        fakeinfile = tmpinfile.replace('.root','')
        tmpoutfile = fakeinfile+year+magtype+"tmp1.root"
        pidvarstokeep = []
        for track, subst in tracks.iteritems() : 
            for var, config in subst.iteritems() : 
                #if building from nightlies
                #command  = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py" 
                command  = "python "+UraniaDir+"/PIDCalib/PIDPerfScripts/scripts/python/PIDGenUser/PIDGen.py"
                command += " -m %s_%s" % (track, ptvar)
                command += " -q %s_%s" % (track, pvar)
                command += " -n %s" % ntrvar
                command += " -t %s" % treename
                command += " -p %s_%s_corr" % (track, var)
                command += " -c %s" % config
                command += " -d %s" % dataset
                command += " -i %s" % tmpinfile
                command += " -o %s" % tmpoutfile
                if seed : 
                  command += " -s %d" % seed
    
                treename  = output_tree
                tmpinfile = tmpoutfile
                if tmpoutfile == fakeinfile+year+magtype+"tmp1.root":
                    tmpoutfile = fakeinfile+year+magtype+"tmp2.root"
                else: 
                    tmpoutfile = fakeinfile+year+magtype+"tmp1.root"
                
                print(command)
                os.system(command)
                pidvarstokeep += ["%s_%s_corr"%(track, var)]
        
        make_tmpfile(tmpinfile, treename, output_file, ['runNumber', 'eventNumber']+pidvarstokeep, nentries_to_read = 1000000000)
        shutil.rmtree(tmpdir)

#def AddPIDCalibWeights(infname, intreename, outfname, magtype, UraniaDir, year = '2016', tmpdir = '/tmp/', perfhistpath = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist/RLc', nentries_to_read = 1000000000):
#    """
#    Function:
#    ---------
#    Adds PIDCalib weights for following cuts
#    Proton: [DLLp>0]
#    Kaon  : [DLLK>4]
#    Pion  : [DLLK<2]
#    Muon  : [DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1] (NB: This is an 'offline' cut on Muon not Stripping cut).
#
#    Arguments:
#    ---------
#    infname:    full path to input file tfile including it's name
#    intreename: ttree name in the input file
#    outfname:   full path to the output file including it's name where it will be written
#    magtype :   MagDown or MagUp
#    year :      2016
#    tmpdir :    Currently default is set to /tmp/.The program makes a tmp directory inside tmpdir i.e. /tmp/tmp. It also later deltes it.
#    UraniaDir:  Path to the Urania directory that holds the run bash script
#    perfhistpath: path where the perfomance histogram are stored (these are created using custom binning scheme). Default is set to '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist/RLc'.
#    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000
#
#    Output:
#    ------
#    Stores a root file with the same TTree name as the input. The file contains eventNumber, runNumber and PIDCalib weights.
#    """
#    Exceptions(infname, year, magtype, tmpdir, UraniaDir)
#
#    if year == '2016': stripp = "Turbo16"
#    print('Using options', infname, intreename, outfname, magtype, year, UraniaDir, perfhistpath)
#
#    tmpdir = make_tmpdir(tmpdir)
#    if not UraniaDir.endswith('/'): UraniaDir += '/'
#    if not perfhistpath.endswith('/'): perfhistpath += '/'
#    pidscript = UraniaDir+'PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/PerformMultiTrackCalib.py'
#
#    tmpfinfname = tmpdir+(infname.split('/')[-1]).replace('.root', '_pidcalibtmp.root')
#    tmpfoutfname = tmpfinfname.replace('.root', '_pidcalibtmp_withWeights.root')
#    leafskeep   = [trck+'_'+lf for trck in trcks for lf in leafs[:2]] + ['runNumber', 'eventNumber', leafs[2]]
#    make_tmpfile(infname, intreename, tmpfinfname, leafskeep, nentries_to_read)
#
#    command = r'bash '+UraniaDir+'run python '+pidscript+' '+stripp+' '+magtype+' '+tmpfinfname+' '+intreename+' '+tmpfoutfname+' \[mu,Mu,DLLmu\>2\&\&DLLmu-DLLK\>2\&\&DLLmu-DLLp\>2\&\&IsMuon\=\=1\] \[pi,Pi,DLLK\<2\] \[K,K,DLLK\>4\] \[p,P,DLLp\>0\] -i '+perfhistpath+' -X Brunel_P -Y Brunel_PT -Z nTracks_Brunel -x P -y PT -z nTracks -s Mu binning-Mu-'+stripp+'-'+magtype+' -s Pi binning-Pi-'+stripp+'-'+magtype+' -s K binning-K-'+stripp+'-'+magtype+' -s P binning-P-'+stripp+'-'+magtype+' -q' 
#    os.system(command)
#
#    make_tmpfile(tmpfoutfname, intreename, outfname, ['runNumber', 'eventNumber', 'Event_PIDCalibEffWeight'], nentries_to_read = 1000000000)
#    shutil.rmtree(tmpdir)

def storeeff(row, histg, limits):
    xmin, xmax, ymin, ymax, zmin, zmax = limits
    glbbin = histg.FindBin(row[0],row[1],row[2]) #global bin number
    cond   = (row[0] < xmin or row[0] > xmax or row[1] < ymin or row[1] > ymax or row[2] < zmin or row[2] > zmax)
    eff    = histg.GetBinContent(glbbin)
    if cond: return 0.

    #if eff < 0.: 
    #    return 0.
    #elif eff > 8.:
    #    return 1.

    return eff

def AddPIDCalibWeights(infname, intreename, outfname, magtype, year = '2016', perfhistpath = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist/RLc', nentries_to_read = 1000000000, chunksize = 10000):
    """
    Function:
    ---------
    Adds PIDCalib weights for following cuts
    Proton: [DLLp>0]
    Kaon  : [DLLK>4]
    Pion  : [DLLK<2]
    Muon  : [DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1&&MC15TuneV1_ProbNNghost<0.2] (NB: This is an 'offline' cut on Muon not Stripping cut).

    Arguments:
    ---------
    infname:    full path to input file tfile including it's name
    intreename: ttree name in the input file
    outfname:   full path to the output file including it's name where it will be written
    magtype :   MagDown or MagUp
    year :      2016
    perfhistpath: path where the perfomance histogram are stored (these are created using custom binning scheme). Default is set to '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist/RLc'.
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000
    chunksize: Pandas data frame chunksize to read

    Output:
    ------
    Stores a root file with the same TTree name as the input. The file contains eventNumber, runNumber and PIDCalib weights (Event_PIDCalibEffWeight).
    """
    TH1.AddDirectory(kFALSE)
    yr        = "Turbo"+year[-2:]
    varsdf    = ['runNumber', 'eventNumber', 'nTracks']
    perfHist  = {}
    for trck_PIDCalib in trcks_PIDCalib:
        if trck_PIDCalib == 'K':
            varsdf     += [trck_PIDCalib+"_P", trck_PIDCalib+"_PT"]
        else:
            varsdf     += [trck_PIDCalib.lower()+"_P", trck_PIDCalib.lower()+"_PT"]

        prefix      = perfhistpath+"/PerfHists_"+trck_PIDCalib+"_"+yr+"_"+magtype
        binningname = "binning-"+trck_PIDCalib+"-"+yr+"-"+magtype
        suffix      = "_".join(leafs_PIDCalib)
        perfname    = prefix+"_"+binningname+"_"+suffix+".root"
        File        = TFile.Open(perfname, "read")
        Histg       = File.Get(perfhistname[trck_PIDCalib])
        perfHist[trck_PIDCalib]  = Histg.Clone(trck_PIDCalib+"new")
        File.Close()

    varstoStore = ['runNumber', 'eventNumber', 'Event_PIDCalibEffWeight']
    if os.path.exists(outfname): os.remove(outfname)
    if nentries_to_read <= chunksize: chunksize = nentries_to_read 
    events_read = 0
    for df_data in read_root(infname, intreename, chunksize=chunksize, columns=varsdf):
        print('Events read', events_read)
        if events_read >= nentries_to_read: break
        for trck_PIDCalib in trcks_PIDCalib:
            Xmin  = perfHist[trck_PIDCalib].GetXaxis().GetXmin(); Xmax  = perfHist[trck_PIDCalib].GetXaxis().GetXmax()
            Ymin  = perfHist[trck_PIDCalib].GetYaxis().GetXmin(); Ymax  = perfHist[trck_PIDCalib].GetYaxis().GetXmax()
            Zmin  = perfHist[trck_PIDCalib].GetZaxis().GetXmin(); Zmax  = perfHist[trck_PIDCalib].GetZaxis().GetXmax()
            Limits= (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)
            if trck_PIDCalib == 'K':
                applyvars = [trck_PIDCalib+'_P', trck_PIDCalib+'_PT', 'nTracks']
                df_data[trck_PIDCalib+'_PIDCalibeff'] = df_data[applyvars].apply(storeeff, args=[perfHist[trck_PIDCalib], Limits], axis=1)
            else:
                applyvars = [trck_PIDCalib.lower()+'_P', trck_PIDCalib.lower()+'_PT', 'nTracks']
                df_data[trck_PIDCalib.lower()+'_PIDCalibeff'] = df_data[applyvars].apply(storeeff, args=[perfHist[trck_PIDCalib], Limits], axis=1)

        df_data['Event_PIDCalibEffWeight'] = df_data['K_PIDCalibeff'] * df_data['p_PIDCalibeff'] * df_data['pi_PIDCalibeff'] * df_data['mu_PIDCalibeff']
        df_data[varstoStore].to_root(outfname, key = intreename, mode='a', store_index=False)
        events_read += df_data.shape[0]
