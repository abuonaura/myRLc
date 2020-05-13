###################
import sys, os
import ROOT as r
import pandas as pd
import numpy as np
from root_pandas import read_root
from sklearn.model_selection import train_test_split
from datetime import datetime
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report, roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import xgboost as xgb
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from matplotlib import gridspec
import joblib
import pickle
###################

###################
datadir   = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/'
files     = {
              'bkg':{'MagUp':datadir+'Data/Data_MagUp_2016.root','MagDown':datadir+'Data/Data_MagDown_2016.root'},
              'sig':datadir+'MC/BDTsamples/Lb2Lctaunu_2016_PIDCalibAndPIDGen.root'
            }
BDTfiles ={ 
               'bkg':{'MagUp': files['bkg']['MagUp'][0:-5]+'_forBDT.root', 'MagDown': files['bkg']['MagDown'][0:-5]+'_forBDT.root'},
               'sig'  : files['sig'][0:-5]+'_forBDT.root'
          }
varsData = ['Lc_M','Lc_FDCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','Lc_IPCHI2_OWNPV','p_ProbNNp','pi_ProbNNpi','K_ProbNNk', 'Lb_L0Global_TIS', 'Lc_L0HadronDecision_TOS', 'Lc_Hlt1TrackMVADecision_TOS', 'Lc_Hlt1TwoTrackMVADecision_TOS', 'Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS', 'mu_PIDmu', 'mu_PIDK', 'mu_PIDp', 'mu_isMuon', 'mu_P', 'pi_P','p_P','K_P', 'mu_PT', 'nTracks']
varsMC   = ['Lc_M', 'Lc_FDCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','Lc_IPCHI2_OWNPV','p_PROBNNp_corr','pi_PROBNNpi_corr','K_PROBNNK_corr', 'Lc_BKGCAT', 'Lb_L0Global_TIS', 'Lc_L0HadronDecision_TOS', 'Lc_Hlt1TrackMVADecision_TOS', 'Lc_Hlt1TwoTrackMVADecision_TOS', 'Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS', 'Lb_TRUEID', 'Lc_TRUEID', 'p_TRUEID', 'K_TRUEID', 'pi_TRUEID', 'mu_TRUEID', 'Event_PIDCalibEffWeight']
BDTvars  = ['Lc_FDCHI2_OWNPV','Lc_IPCHI2_OWNPV','Lc_ENDVERTEX_CHI2','Sum_Lcdaug_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','p_ProbNNp','pi_ProbNNpi','K_ProbNNk']
###################

###################
def createBDTsamples():
    #Apply some preselection: Trigger + Lc_M + for MC Truth matching. Also look at ApplyCuts function
    cut    = "(Lb_L0Global_TIS==1||Lc_L0HadronDecision_TOS==1)&&(Lc_Hlt1TrackMVADecision_TOS==1||Lc_Hlt1TwoTrackMVADecision_TOS==1)&&Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS==1&&Lc_M>2230&&Lc_M<2330"
    cut_mc = cut+"&&abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122&&abs(p_TRUEID)==2212&&abs(K_TRUEID)==321&&abs(pi_TRUEID)==211&&abs(mu_TRUEID)==13"

    #Make a copy of MC for BDT traning
    f_sig = r.TFile(files['sig'],'READ')
    t_sig = f_sig.Get('DecayTree')
    t_sig.SetBranchStatus('*',0)
    for var in varsMC: t_sig.SetBranchStatus(var,1)
    of_sig = r.TFile(BDTfiles['sig'],'recreate')
    ot_sig = r.TTree('DecayTree','DecayTree')
    print('>>> Copying signal file: '+files['sig'])
    #Should inprinciple apply '&&Lc_M>2260&&Lc_M<2310&Lb_BKGCAT<30.' but done later (need full sample for performance studies)
    ot_sig =t_sig.CopyTree(cut_mc)
    print('>>> Created signal file for BDT: '+BDTfiles['sig'])
    ot_sig.Write()
    of_sig.Close()
    of_sig.Close()

    #Make a copy of Data for BDT traning - MagUp
    f_bkg_MagUp = r.TFile(files['bkg']['MagUp'],'READ')
    t_bkg_MagUp = f_bkg_MagUp.Get('tupleout/DecayTree')
    t_bkg_MagUp.SetBranchStatus('*',0)
    for var in varsData: t_bkg_MagUp.SetBranchStatus(var,1)
    ofname_bkg_MagUp   =files['bkg']['MagUp'][0:-5]+'_forBDT.root'
    of_bkg_MagUp = r.TFile(BDTfiles['bkg']['MagUp'],'recreate')
    ot_bkg_MagUp = r.TTree('DecayTree','DecayTree')
    print('>>> Copying bkg file: '+files['bkg']['MagUp'])
    #Should inprinciple apply '&&(Lc_M<2260. || Lc_M>2310.)' but done later (need full sample for performance studies)
    ot_bkg_MagUp =t_bkg_MagUp.CopyTree(cut,'',1000000) #only read 100000 in sideband
    print('>>> Created bkg file for BDT: '+BDTfiles['bkg']['MagUp'])
    ot_bkg_MagUp.Write()
    of_bkg_MagUp.Close()
    f_bkg_MagUp.Close()

    #Make a copy of Data for BDT traning - MagDown
    f_bkg_MagDown = r.TFile(files['bkg']['MagDown'],'READ')
    t_bkg_MagDown = f_bkg_MagDown.Get('tupleout/DecayTree')
    t_bkg_MagDown.SetBranchStatus('*',0)
    for var in varsData:t_bkg_MagDown.SetBranchStatus(var,1)
    of_bkg_MagDown = r.TFile(BDTfiles['bkg']['MagDown'],'recreate')
    ot_bkg_MagDown = r.TTree('DecayTree','DecayTree')
    print('>>> Copying bkg file: '+files['bkg']['MagDown'])
    #Should inprinciple apply '&&(Lc_M<2260. || Lc_M>2310.)' but done later (need full sample for performance studies)
    ot_bkg_MagDown =t_bkg_MagDown.CopyTree(cut,'',1000000) #only read 100000 in sideband
    print('>>> Created bkg file for BDT: '+BDTfiles['bkg']['MagDown'])
    ot_bkg_MagDown.Write()
    of_bkg_MagDown.Close()
    f_bkg_MagDown.Close()

    return None

def ApplyCuts(df):
    cond   = ((df['mu_PIDmu'] > 2) & ((df['mu_PIDmu']-df['mu_PIDK']) > 2) & ((df['mu_PIDmu']-df['mu_PIDp']) > 2) & (df['mu_isMuon']==1))
    cond   = cond & (df['pi_PT'] > 0) & (df['pi_PT'] < 60000) & (df['pi_P'] > 0) & (df['pi_P'] < 200000) 
    cond   = cond & (df['K_PT']  > 0) & (df['K_PT']  < 60000) & (df['K_P']  > 0) & (df['K_P']  < 200000) 
    cond   = cond & (df['p_PT']  > 0) & (df['p_PT']  < 60000) & (df['p_P']  > 0) & (df['p_P']  < 200000) 
    cond   = cond & (df['mu_PT'] > 0) & (df['mu_PT'] < 60000) & (df['mu_P'] > 0) & (df['mu_P'] < 200000) 
    cond   = cond & (df['nTracks'] > 0)&(df['nTracks'] < 700)
    return cond

def CreateTrainTest():
    #Get the files - 
    df_data_up   = read_root(BDTfiles['bkg']['MagUp'],  'DecayTree', columns=varsData)
    df_data_up["Event_PIDCalibEffWeight"] = np.array(df_data_up.shape[0]*[1.])
    df_data_up["Sum_Lcdaug_PT"] = df_data_up["K_PT"] + df_data_up["pi_PT"] + df_data_up["p_PT"]
    #(PIDCalibw include cut on muon vars, other variables are already cut on in the stripping)
    cond       = ApplyCuts(df_data_up)
    df_data_up = df_data_up[cond]
    df_data_down = read_root(BDTfiles['bkg']['MagDown'],'DecayTree', columns=varsData)
    df_data_down["Event_PIDCalibEffWeight"] = np.array(df_data_down.shape[0]*[1.])
    df_data_down["Sum_Lcdaug_PT"] = df_data_down["K_PT"] + df_data_down["pi_PT"] + df_data_down["p_PT"]
    cond         = ApplyCuts(df_data_down)
    df_data_down = df_data_down[cond]

    df_sig       = read_root(BDTfiles['sig'],'DecayTree', columns=varsMC)
    df_sig = df_sig.rename(columns={"p_PROBNNp_corr": "p_ProbNNp","pi_PROBNNpi_corr": "pi_ProbNNpi","K_PROBNNK_corr" : "K_ProbNNk"})
    df_sig["Sum_Lcdaug_PT"] = df_sig["K_PT"] + df_sig["pi_PT"] + df_sig["p_PT"]

    #Select as background only events in Lc sidebands and put together he 2 samples
    df_bkg_up   = df_data_up
    df_bkg_down = df_data_down
    print(df_bkg_up.shape)
    print(df_bkg_down.shape)
    print(df_sig.shape)
    df_bkg_up   = df_bkg_up.loc[(df_bkg_up['Lc_M']<2260)|(df_bkg_up['Lc_M']>2310.)].dropna()
    df_bkg_down = df_bkg_down.loc[(df_bkg_down['Lc_M']<2260)|(df_bkg_down['Lc_M']>2310.)].dropna()
    df_bkg = pd.concat([df_bkg_up, df_bkg_down])
    #Select as signal only MC events in 3sigma eround the Lc mass peak and BKGCAT<30
    df_sig = df_sig.loc[(df_sig['Lc_M']>2260)&(df_sig['Lc_M']<2310.)&(df_sig['Lc_BKGCAT']<30.)].dropna()
    print(df_bkg_up.shape)
    print(df_bkg_down.shape)
    print(df_sig.shape)

    #Do the train and test split sample
    df_bkg_red = df_bkg[:100000] #Out of 100000 used only 20000 events in the training
    X = np.concatenate((df_sig[BDTvars+["Event_PIDCalibEffWeight"]], df_bkg_red[BDTvars+["Event_PIDCalibEffWeight"]]))
    y = np.concatenate((np.ones(df_sig.shape[0]),np.zeros(df_bkg_red.shape[0]))) #Sigs: 1 and Bkgs: 0
    y = y.reshape(y.shape[0], -1)
    df = pd.DataFrame(np.hstack((X,y)),columns=BDTvars+["Event_PIDCalibEffWeight"]+['y'])
    dfX = df[BDTvars+["Event_PIDCalibEffWeight"]]
    dfY = df['y']
    X_train, X_test, y_train, y_test = train_test_split(dfX, dfY, test_size=0.4, random_state=0)
    print('Number of events in the training sample:  ', X_train.shape[0])
    print('Number of events in the test sample:  ', X_test.shape[0])
    ws_train = X_train['Event_PIDCalibEffWeight']
    ws_test  = X_test['Event_PIDCalibEffWeight']
    ws_train_sig = X_train['Event_PIDCalibEffWeight'][y_train>0.5]
    ws_train_bkg = X_train['Event_PIDCalibEffWeight'][y_train<0.5]
    ws_test_sig  = X_test['Event_PIDCalibEffWeight'][y_test>0.5]
    ws_test_bkg  = X_test['Event_PIDCalibEffWeight'][y_test<0.5]
    X_train = X_train[BDTvars]
    X_test  = X_test[BDTvars]
    ws = (ws_train, ws_test, ws_train_sig, ws_train_bkg, ws_test_sig, ws_test_bkg)

    ############ - Plot input vars and LcM
    #Plot the Lc mass distribution for signal+bkg
    fig = plt.figure(1)
    bins = np.linspace(2230, 2330, 50)
    plt.hist(df_bkg['Lc_M'], bins, density=1, facecolor='red', alpha=0.75, label = 'Bkg evts for BDT',edgecolor='k', weights=df_bkg['Event_PIDCalibEffWeight'])
    plt.hist(df_sig['Lc_M'], bins, density=1, facecolor='blue', alpha=0.75, label = 'Signal evts for BDT',edgecolor='k', weights=df_sig['Event_PIDCalibEffWeight'])
    plt.xlabel('$\Lambda_{c}$ mass (MeV)')
    plt.legend(prop={'size':10})
    #plt.show(block = False)
    plt.savefig("plots/LcM_BDTfig.pdf")

    #Decide to check the distributions of the different BDT variables for signal/bkg datasets
    fig1 = plt.figure(2)
    gs = gridspec.GridSpec(3, 4)
    for i,var in enumerate(BDTvars):
            ax = fig1.add_subplot(gs[i])
            low = min(df_sig[var].min(), df_bkg[var].min())
            high = max(df_sig[var].max(), df_bkg[var].max())
            ax.hist(df_sig[var].dropna().values, density =1,color = 'blue',range=(low,high), alpha=0.5, weights=df_sig['Event_PIDCalibEffWeight'])
            ax.hist(df_bkg[var].dropna().values, density=1,color = 'red',range=(low,high), alpha=0.5, weights=df_bkg['Event_PIDCalibEffWeight'])
            ax.set_yscale("log", nonposy='clip')
            if var=='Lc_FDCHI2_OWNPV' or var=='Lc_IPCHI2_OWNPV':
                ax.set_xscale("log")
            ax.set_title(var)
    plt.subplots_adjust(hspace = 1., wspace=1.)
    #plt.show(block=False)
    plt.savefig("plots/Inpurvar_BDTfig.pdf")
    ############
    return df_data_up, df_data_down, df_sig, X_train, X_test, y_train, y_test, ws

def FitModel(X_train, X_test, y_train, y_test, sample_wghts=None):
    #define fit model with optimised hyperparameters
    model = XGBClassifier(base_score=0.5, booster='gbtree', colsample_bylevel=1,
           colsample_bytree=1.0, gamma=2, learning_rate=0.02,
           max_delta_step=0, max_depth=5, min_child_weight=0.1, missing=None,
           n_estimators=600, n_jobs=1, nthread=1, objective='binary:logistic',
           random_state=0, reg_alpha=0, reg_lambda=1, scale_pos_weight=1,
           seed=None, silent=True, subsample=0.6)

    #fit model
    model.fit(X_train, y_train, sample_weight=sample_wghts)

    #########
    # Performance on train
    y_pred_train = model.predict(X_train)
    auc_score = roc_auc_score(y_train, y_pred_train)
    print("Performance on train sample roc_auc_score: ", auc_score)

    # Performance on test
    y_pred_test = model.predict(X_test)
    predictions = [round(value) for value in y_pred_test]
    auc_score = roc_auc_score(y_test, y_pred_test)
    accuracy = accuracy_score(y_test, predictions)
    print("Performance on test sample roc_auc_score: ", auc_score)
    print("accuracy_score: %.2f%%" % (accuracy * 100.0))
    #########

    return model, y_pred_test

def RunBDT(BDTcut):
    if not os.path.isdir('plots'): os.mkdir('plots')
    createBDTsamples()
    df_data_up, df_data_down, df_sig, X_train, X_test, y_train, y_test, ws = CreateTrainTest()
    model, y_pred_test = FitModel(X_train, X_test, y_train, y_test, sample_wghts=ws[0])
    print('Using predict_proba')
    fpr, tpr, areaUC = PrintROC(X_test, y_test, y_pred_test, model, 'xgboost')
    decisions        = compare_train_test(model, X_train, y_train, X_test, y_test)
    plot_compare_train_test(decisions,25,'xgboost', ws = ws[2:])
    PlotSingleROC(fpr, tpr, areaUC, 'xgboost', 'navy')
    PlotFeatureImportance(model,'xgboost')
    #To see what are the effects of BDT on Lcmass let's apply the trained model onto the full data dataset
    df_data = pd.concat([df_data_up, df_data_down])
    y_predicted = model.predict_proba(df_data[BDTvars])[:,1]
    df_data     = df_data.assign(BDT=y_predicted)
    PlotLcMass(df_data, BDTcut)
    file_name = "xgb_reg.pkl"
    pickle.dump(model, open("xgb_reg.pkl", "wb"))

    return None

def LoadPIDGenDataF(infname):
    df_list = []
    fname = infname[0:-5]+'_PIDGen.root'
    for df in read_root(fname, 'DecayTree', chunksize=10000):
        df_list.append(df)
    return df_list

def AddBDTinfo(infname, intreename, outfname, datatype, pickled_model_path = './xgb_reg.pkl', nentries_to_read=1000000000):
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
    pickled_model_path: full path to the saved tranined model (it must be pickled). Default is './xgb_reg.pkl'
    nentries_to_read: number of candidates to read. Set to root default i.e. 1000000000

    Output:
    ------
    Stores a root file with the same TTree name as the input. The file contains eventNumber, runNumber and PIDCalib weights.
    """
    model = pickle.load(open(pickled_model_path, "rb"))
    outtreename = intreename.split('/')[-1]
    if datatype == 'MC': 
        df_PIDgen = LoadPIDGenDataF(infname)
        i=0
    for df_data in read_root(infname, intreename, chunksize=10000):
        df_forpred = None
        if datatype == 'MC': 
            df_merge = df_data.merge(df_PIDgen[i], left_index=True, right_index=True)
            i+=1
            #print(varsMC[:13])
            df_forpred = df_merge[varsMC[:13]]
            df_forpred = df_forpred.rename(columns={"p_PROBNNp_corr": "p_ProbNNp","pi_PROBNNpi_corr": "pi_ProbNNpi","K_PROBNNK_corr" : "K_ProbNNk"})
        else:
            df_forpred = df_data[varsData]

        df_forpred["Sum_Lcdaug_PT"] = df_forpred["K_PT"] + df_forpred["pi_PT"] + df_forpred["p_PT"]
        y_predicted = model.predict_proba(df_forpred[BDTvars])[:, 1]
        df_data['bdt'] = y_predicted
        #for d_k in list(df_data.keys()): df_data[d_k] = pd.to_numeric(df_data[d_k], errors='coerce')
        #print(df_data.shape)
        df_data = df_data[['eventNumber', 'runNumber', 'bdt']]
        df_data.to_root(outfname, key = outtreename, mode='a', store_index=False)
        #print(df_data.shape[0])
        if df_data.shape[0] > nentries_to_read: 
            print('Asked to only read', nentries_to_read, ', however I have read', df_data.shape[0])
            break

    #f = r.TFile(fname,'update')
    #if not f: print('File NOT found!')
    #t = f.Get('DecayTree')
    #if not t: print('Tree NOT found!')
    #bdt       = np.zeros(1, dtype=float)
    #bdtBranch = t.Branch("bdt",bdt,"bdt/D")
    #if dtype!='MC':
    #    df_data = read_root(fname, 'DecayTree', columns=varsData)
    #else:
    #    df_data = read_root(fname, 'DecayTree', columns=varsMC)
    #    df_data = df_data.rename(columns={"p_PROBNNp_corr": "p_ProbNNp","pi_PROBNNpi_corr": "pi_ProbNNpi","K_PROBNNK_corr" : "K_ProbNNk"})
    #y_predicted = model.predict_proba(df_data[BDTvars])[:, 1]
    #for i in range(t.GetEntries()):
    #    t.GetEntry(i)
    #    bdt[0] = y_predicted[i]
    #    bdtBranch.Fill()
    #f.Write()
    #f.Close()
    return None

def timer(start_time=None):
    if not start_time:
        start_time = datetime.now()
        return start_time
    elif start_time:
        thour, temp_sec = divmod((datetime.now() - start_time).total_seconds(), 3600)
        tmin, tsec = divmod(temp_sec, 60)
        print('\n Time taken: %i hours %i minutes and %s seconds.' % (thour, tmin, round(tsec, 2)))

def PerformGridSearch(X_train, y_train):
    # A parameter grid for XGBoost
    params = {
                'min_child_weight': [0.1, 0.5, 1, 5],
                'gamma': [0.5, 1, 1.5, 2],
                'subsample': [0.3, 0.6, 1.0],
                'colsample_bytree': [0.6, 0.8, 1.0],
                'max_depth': [3, 4, 5]
                }
    xgb = XGBClassifier(learning_rate=0.02, n_estimators=600,silent=True, nthread=1)
    folds = 5
    param_comb = 100
    skf = StratifiedKFold(n_splits=folds, shuffle = True, random_state = 1001)
    random_search = RandomizedSearchCV(xgb, param_distributions=params, n_iter=param_comb, scoring='roc_auc', n_jobs=10, cv=skf.split(X_train,y_train), verbose=3, random_state=1001 )
    start_time = timer(None) # timing starts from this point for "start_time" variable
    random_search.fit(X_train, y_train)
    timer(start_time) # timing ends here for "start_time" variable
    print('\n All results:')
    print(random_search.cv_results_)
    print('\n Best estimator:')
    print(random_search.best_estimator_)
    print('\n Best normalized gini score for %d-fold search with %d parameter combinations:' % (folds, param_comb))
    print(random_search.best_score_ * 2 - 1)
    print('\n Best hyperparameters:')
    print(random_search.best_params_)
    results = pd.DataFrame(random_search.cv_results_)
    results.to_csv('xgb-random-grid-search-results-01.csv', index=False)
    return None

def PrintROC(X_test, y_test, y_predicted, model, classifier_name):
    fpr_lr, tpr_lr, thr = roc_curve(y_test, model.predict_proba(X_test)[:, 1])
    roc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
    areaUC = auc(fpr_lr, tpr_lr)
    print("Classifier used: ", classifier_name)
    print(classification_report(y_test, y_predicted, target_names=["background", "signal"]))
    print("Area under ROC curve: %.4f"%(roc))
    return fpr_lr, tpr_lr, areaUC

def compare_train_test(model, X_train, y_train, X_test, y_test) :
    """
    Given a trained classifier, get the distribution of the classifier output for the signal and background categories of the training and testing datasets.
    This allows to check for overtraining.
    The output of this function can be passed to other functions for python (python_plot_compare_train_test)
    """
    decisions = []
    for X,y in ((X_train, y_train), (X_test, y_test)) :
        predS = model.predict_proba(X[y>0.5])[:, 1].ravel()
        predB = model.predict_proba(X[y<0.5])[:, 1].ravel()
        decisions += [predS, predB]
    return decisions

def plot_compare_train_test(decisions,bins,classifier, ws=None):
    """
    Plot the distribution of the classifier output for the signal and background categories of the training and testing datasets.
    This allows to check for overtraining.
    """
    low  = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    low_high = (low,high)
    # Plot with python.
    plt.figure()
    plt.hist(decisions[0], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='S (train)', weights=ws[0])
    plt.hist(decisions[1], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='B (train)', weights=ws[1])
    hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True, weights=ws[2])
    center = (bins[:-1] + bins[1:]) / 2
    #scale = len(decisions[2]) / sum(hist)
    scale = sum(ws[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')
    hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True, weights=ws[3])
    #scale = len(decisions[3]) / sum(hist)
    scale = sum(ws[3]) / sum(hist)
    err = np.sqrt(hist * scale) / scale
    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')
    plt.xticks(np.arange(0, 1, step=0.1))
    plt.xlabel("Classifier output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    plt.savefig('plots/plt_' + classifier+'_Output.pdf',format='pdf')
    plt.show(block = False)
    return None

def PlotSingleROC(fpr, tpr, roc_auc, classifier, color):
    #Plot ROC curve
    plt.figure()
    plt.xlim([-0.01, 1.00])
    plt.ylim([-0.01, 1.01])
    plt.plot(fpr, tpr, lw=1, color = color, label='ROC curve (area = %0.3f)' % roc_auc)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC curve '+ classifier, fontsize=14)
    plt.legend(loc='lower right', fontsize=12)
    plt.savefig('plots/ROCcurve_'+classifier+'.png')
    plt.show(block = False)
    return None

def PlotFeatureImportance(model,classifier):
    features=BDTvars
    importances = model.feature_importances_
    indices = np.argsort(importances)
    fig = plt.figure()
    plt.title('Feature Importances')
    ax = fig.add_subplot(111)
    plt.bar(range(len(importances)), importances,color='steelblue')
    ticks = np.arange(0,12,1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(features, rotation = 'vertical')
    plt.savefig('plots/figureImportances'+classifier+'.png',bbox_inches='tight')
    plt.show(block=False)
    return None

def PlotLcMass(df_data, BDTcut, ws=None):
    mLc_s = df_data['Lc_M'].where(df_data['BDT']>BDTcut).dropna()
    mLc_b = df_data['Lc_M'].where(df_data['BDT']<BDTcut).dropna()
    mLc_tot = df_data['Lc_M']
    w_s   = df_data['Event_PIDCalibEffWeight'].where(df_data['BDT']>BDTcut).dropna()
    w_b   = df_data['Event_PIDCalibEffWeight'].where(df_data['BDT']<BDTcut).dropna()
    w_tot = df_data['Event_PIDCalibEffWeight']
    #Plot
    fig2 = plt.figure()
    plt_s = plt.subplot(1,3,1)
    plt_s.set_xlabel('$\Lambda_{c}$ mass (MeV)')
    plt.hist(mLc_s, color = 'blue', bins=30, edgecolor='k')
    plt_s.set_title('BDT > '+ str(BDTcut))
    plt_b = plt.subplot(1,3,2,sharex = plt_s)
    plt_b.set_xlabel('$\Lambda_{c}$ mass (MeV)')
    plt.hist(mLc_b,color='red',bins=30,edgecolor='k')
    plt_b.set_title('BDT < '+ str(BDTcut))
    plt_tot = plt.subplot(1,3,3,sharex = plt_s)
    plt_tot.set_xlabel('$\Lambda_{c}$ mass (MeV)')
    plt_tot.hist(mLc_tot, bins=30, density=False, alpha = 0.5,color = 'green', label='No cut',edgecolor='k', weights=w_tot)
    plt_tot.hist(mLc_b, bins=30, density=False,alpha = 0.5,color='red', label='BDT < '+str(BDTcut),edgecolor='k', weights=w_b)
    plt_tot.hist(mLc_s, bins=30, density=False,alpha = 0.5,color='blue', label='BDT > '+str(BDTcut),edgecolor='k', weights=w_s)
    plt_tot.legend(prop={'size': 10})
    plt.savefig('plots/mLcWithBDT.png',bbox_inches='tight')
    plt.show(block = False)
    return None
###################

#RunBDT(0.7)
