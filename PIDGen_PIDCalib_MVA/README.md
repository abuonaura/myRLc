Scripts for storing PIDGen vars, PIDCalib and MVA info. 
Each stage is factored into a function that stores a root file that contains eventNumber, runNumber and the information regarding that step.
See Store_PIDGen_PIDCalib_MVA_Info.py for the example.

Before running, I recommend copying the following scripts to the locations mentioned.
The rest of the document explains the functions and scripts used in the stages.

# Altered PIDCalib scripts
I have slightly altered the PIDCalib scripts to fit my needs, which have been stored in directory PIDCalib_AlteredScript.
These files need to be copied to your corresponding Urania directory and Urania needs to be built again.

```bash
    cp PIDCalib_AlteredScript/binningPID.py $UraniaDir/PIDCalib/PIDPerfScripts/scripts/python/BinningOptimizer/
```

This script is used for PIDCalib binning optimisation. 
I have modified this script to read only set number of PIDCalib files, default is that it runs over all of them and take quite some time.

```bash
    cp PIDCalib_AlteredScript/PerformMultiTrackCalib.py $UraniaDir/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/
    cp PIDCalib_AlteredScript/MultiTrackCalibTool.cpp $UraniaDir/PIDCalib/PIDPerfTools/src/
```

Modified these scripts to use the ttree branch and the associated leaves of the reference MC sample directly instead of writing a separate root file.
Maybe AddFriend method takes care of this already and these modifications may not be needed.

```bash
    cd $UraniaDir
    make purge; make clean
    make configure; make; make install
```

# Function to store PIDGen vars
Store PIDGen vars. The function resides in file Add_PIDGenWghts_PIDCalibWghts.py and requires Urania installation.

```python
    AddPIDGenWeights(infname, intreename, outfname, magtype, UraniaDir, year = '2016', tmpdir = '/tmp/', nentries_to_read = 1000000000)
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
    tmpdir :    Currently default is set to /tmp/.The program makes a tmp directory inside tmpdir i.e. /tmp/tmp. It also later deletes it.
    UraniaDir:  Path to the Urania directory that holds the run bash script
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000
    
    Output:
    ------
    Stores a root file with the same TTree name as the input. It only contains eventNumber, runNumber and PIDGen vars as branches.
    """
```

# Function to store PIDCalib weights
Store PIDCalib weights. Function is defined in Add_PIDGenWghts_PIDCalibWghts.py and requires performance histograms and built Urania.

```python
    AddPIDCalibWeights(infname, intreename, outfname, magtype, UraniaDir, year = '2016', tmpdir = '/tmp/', perfhistpath = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist', nentries_to_read = 1000000000)
    """
    Function:
    ---------
    Adds PIDCalib weights for following cuts
    Proton: [DLLp>0]
    Kaon  : [DLLK>4]
    Pion  : [DLLK<2]
    Muon  : [DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1] (NB: This is an 'offline' cut on Muon not Stripping cut).

    Arguments:
    ---------
    infname:    full path to input file tfile including it's name
    intreename: ttree name in the input file
    outfname:   full path to the output file including it's name where it will be written
    magtype :   MagDown or MagUp
    year :      2016
    tmpdir :    Currently default is set to /tmp/.The program makes a tmp directory inside tmpdir i.e. /tmp/tmp. It also later deletes it.
    UraniaDir:  Path to the Urania directory that holds the run bash script
    perfhistpath: path where the performance histogram are stored (these are created using custom binning scheme).
    perfhistpath: path where the performance histogram are stored (these are created using custom binning scheme). Default is set to '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist'.
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000

    Output:
    ------
    Stores a root file with the same TTree name as the input. The file contains eventNumber, runNumber and PIDCalib weights as branches.
    """
```

# Function to store BDT info for MC and Data
Store BDT info for MC and Data. Function defined in Add_MVA.py (Adopted from Annarita) and requires a LOT of python packages (See the end on how to install them).
This script can also be used for training MVA but need to set the signal and bkg sample paths correctly.

```python
    AddBDTinfo(infname, intreename, outfname, datatype, pickled_model_path = './xgb_reg.pkl', nentries_to_read=1000000000):
    """
    Function:
    ---------
    Loads the trained pickled MVA model and evaluates on both data and MC samples

    Arguments:
    ---------
    infname:    full path to input file tfile including it's name
    intreename: ttree name in the input file
    outfname:   full path to the output file including it's name where it will be written
    datatype:   MC or Data
    pickled_model_path: full path to the saved trained model (it must be pickled). Default is './xgb_reg.pkl'
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000

    Output:
    ------
    Stores a root file with the same TTree name as the input. The file contains eventNumber, runNumber and PIDCalib weights.
    """
```

# Script to make performance histogram
Make performance histogram for the cut we applied i.e. 
Proton: [DLLp>0]
Kaon  : [DLLK>4]
Pion  : [DLLK<2]
Muon  : [DLLmu>2&&DLLmu-DLLK>2&&DLLmu-DLLp>2&&IsMuon==1] (NB: This is an 'offline' cut on Muon not Stripping cut).

```bash
    bash Make_PerfHist.sh <Urania_path> <path_to_store_perfhist>
    e.g. bash Make_PerfHist.sh /home/hep/amathad/Packages/UraniaDev_v8r0 /disk/lhcb_data2/amathad/Lb2Lclnu_analysis/perfhist
```
Outputs a .root perfhist in the specified directory.

# Script to 1D optimal PIDCalib binning scheme
How to get an optimal 1D PIDCalib binning scheme? 
The following script makes configuration files in direc 'config_files' for each particle, year, magtype, kinematic variable.
Then passes each of the config file to binningPID.py script (see above), which then spits out ONE file (e.g. binning-Turbo16.py) that contains the binning scheme for all the particles, years, magtypes, kinematic variables.

```python
    python Make_ConfigFiles_BinningScheme.py <Urania_path>
    e.g. python Make_ConfigFiles_BinningScheme.py /home/hep/amathad/Packages/UraniaDev_v8r0
```
Outputs 'binning-Turbo16.py' scripts that contains the optimal binning scheme to be used for all particles, years, magtypes and kinematic variables. 

# Figure of merit optimisation
Script for MVA cut optimisation. Correct path for Data and MC need to be set.

```python
    python FOM_optimisation.py
```
Outputs 'FOMcourse/FOM.csv' containing columns 'MVACutval:Sig#;Bkg#:eff_mva:S/Sqrt(S+B):S^2/(S+B)^(3/2)'

# Installing packages to run xGBoost
How to install all the packages to run xGBoost?

```bash
    conda create -n tfrootenv python=3.7 root -c conda-forge
    conda activate tfrootenv
    #conda config --env --add channels conda-forge
    #conda config --add channels conda-forge
    conda install root_numpy root_pandas scikit-learn ipython pandas matplotlib scipy scikit-learn sympy seaborn tensorflow-gpu=1.13.1
    conda install -c conda-forge xgboost tensorflow-probability==0.6.0
    conda install -c conda-forge zfit==0.3.5
```
