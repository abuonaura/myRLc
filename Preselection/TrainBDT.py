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

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
files = {'bkg':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced.root','MagDown':datadir+'Data/Lb_Data_MagDown_reduced.root'},'sig':datadir+'MC/Lb_Lctaunu_PID_reduced.root'}

varsData = ['Lc_M','Lc_FDCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','Lc_IPCHI2_OWNPV','p_ProbNNp','pi_ProbNNpi','K_ProbNNk']

varsMC = ['Lc_M','Lc_BKGCAT','Lc_FDCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','Lc_IPCHI2_OWNPV','p_ProbNNp_corr','pi_ProbNNpi_corr','K_ProbNNk_corr']

BDTvars =['Lc_FDCHI2_OWNPV','Lc_IPCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','p_ProbNNp','pi_ProbNNpi','K_ProbNNk']

def createBDTsamples():
    f_sig = r.TFile(files['sig'],'READ')
    t_sig = f_sig.Get('DecayTree')
    t_sig.SetBranchStatus('*',0)
    for var in varsMC:
        t_sig.SetBranchStatus(var,1)
    #cut_sig ='Lc_M>2260. && Lc_M<2310. && Lc_BKGCAT<30.'


    ofname_sig =files['sig'][0:-5]+'_forBDT.root'
    of_sig = r.TFile(ofname_sig,'recreate')
    ot_sig = r.TTree('DecayTree','DecayTree')
    print('>>> Copying signal file: '+files['sig'])
    ot_sig =t_sig.CopyTree('')
    print('>>> Created signal file for BDT: '+ofname_sig)
    ot_sig.Write()
    of_sig.Close()
    of_sig.Close()

    f_bkg_MagUp = r.TFile(files['bkg']['MagUp'],'READ')
    t_bkg_MagUp = f_bkg_MagUp.Get('DecayTree')
    t_bkg_MagUp.SetBranchStatus('*',0)
    for var in varsData:
        t_bkg_MagUp.SetBranchStatus(var,1)
    #cut_bkg = 'Lc_M<2260. || Lc_M>2310.'

    ofname_bkg_MagUp =files['bkg']['MagUp'][0:-5]+'_forBDT.root'
    ofname_bkg_MagDown =files['bkg']['MagDown'][0:-5]+'_forBDT.root'
    of_bkg_MagUp = r.TFile(ofname_bkg_MagUp,'recreate')
    ot_bkg_MagUp = r.TTree('DecayTree','DecayTree')
    print('>>> Copying bkg file: '+files['bkg']['MagUp'])
    ot_bkg_MagUp =t_bkg_MagUp.CopyTree('','',100000)
    print('>>> Created bkg file for BDT: '+ofname_bkg_MagUp)
    ot_bkg_MagUp.Write()
    of_bkg_MagUp.Close()
    f_bkg_MagUp.Close()

    f_bkg_MagDown = r.TFile(files['bkg']['MagDown'],'READ')
    t_bkg_MagDown = f_bkg_MagDown.Get('DecayTree')
    t_bkg_MagDown.SetBranchStatus('*',0)
    for var in varsData:
        t_bkg_MagDown.SetBranchStatus(var,1)

    of_bkg_MagDown = r.TFile(ofname_bkg_MagDown,'recreate')
    ot_bkg_MagDown = r.TTree('DecayTree','DecayTree')
    print('>>> Copying bkg file: '+files['bkg']['MagDown'])
    ot_bkg_MagDown =t_bkg_MagDown.CopyTree('','',100000)
    print('>>> Created bkg file for BDT: '+ofname_bkg_MagDown)
    ot_bkg_MagDown.Write()
    of_bkg_MagDown.Close()
    f_bkg_MagDown.Close()
    return


def CreateTrainTest():

    BDTfiles={'Data':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced_forBDT.root',
        'MagDown':datadir+'Data/Lb_Data_MagDown_reduced_forBDT.root'},
        'MC':datadir+'MC/Lb_Lctaunu_PID_reduced_forBDT.root'}

    df_data_up = read_root(BDTfiles['Data']['MagUp'],'DecayTree', columns=varsData)
    df_data_down = read_root(BDTfiles['Data']['MagDown'],'DecayTree', columns=varsData)
    df_sig = read_root(BDTfiles['MC'],'DecayTree', columns=varsMC)

    #Rename some variables so that they have the same name for bkg
    df_sig =df_sig.rename(columns={"p_ProbNNp_corr": "p_ProbNNp","pi_ProbNNpi_corr": "pi_ProbNNpi","K_ProbNNk_corr" : "K_ProbNNk"})

    #Select as background only events in Lc sidebands and put together he 2 samples
    df_bkg_up = df_data_up
    df_bkg_down = df_data_down
    df_bkg_up = df_bkg_up.loc[(df_bkg_up['Lc_M']<2260)|(df_bkg_up['Lc_M']>2310.)].dropna()
    df_bkg_down = df_bkg_down.loc[(df_bkg_down['Lc_M']<2260)|(df_bkg_down['Lc_M']>2310.)].dropna()
    frames = [df_bkg_up, df_bkg_down]
    df_bkg = pd.concat(frames)
    #Select as signal only MC events in 3sigma eround the Lc mass peak and BKGCAT<30
    df_sig = df_sig.loc[(df_sig['Lc_M']>2260)&(df_sig['Lc_M']<2310.)&(df_sig['Lc_BKGCAT']<30.)].dropna()

    #Plot the Lc mass distribution for signal+bkg
    fig = plt.figure(1)
    bins = np.linspace(2230, 2330, 50)
    plt.hist(df_bkg['Lc_M'], bins, density=1, facecolor='red', alpha=0.75, label = 'Bkg evts for BDT',edgecolor='k')
    plt.hist(df_sig['Lc_M'], bins, density=1, facecolor='blue', alpha=0.75, label = 'Signal evts for BDT',edgecolor='k')
    plt.xlabel('$\Lambda_{c}$ mass (MeV)')
    plt.legend(prop={'size':10})
    #plt.show(block = False)

    #Decide to check the distributions of the different BDT variables for signal/bkg datasets
    fig1 = plt.figure(2)
    gs = gridspec.GridSpec(3, 4)
    for i,var in enumerate(BDTvars):
            ax = fig1.add_subplot(gs[i])
            low = min(df_sig[var].min(), df_bkg[var].min())
            high = max(df_sig[var].max(), df_bkg[var].max())
            ax.hist(df_sig[var].dropna().values, density =1,color = 'blue',range=(low,high), alpha=0.5)
            ax.hist(df_bkg[var].dropna().values, density=1,color = 'red',range=(low,high), alpha=0.5)
            ax.set_yscale("log", nonposy='clip')
            if var=='Lc_FDCHI2_OWNPV' or var=='Lc_IPCHI2_OWNPV':
                ax.set_xscale("log")
            ax.set_title(var)
    plt.subplots_adjust(hspace = 1., wspace=1.)
    plt.show(block=False)

    df_bkg_red = df_bkg[:20000]
    X = np.concatenate((df_sig[BDTvars], df_bkg_red[BDTvars]))
    y = np.concatenate((np.ones(df_sig.shape[0]),np.zeros(df_bkg_red.shape[0])))
    y = y.reshape(y.shape[0], -1)
    df = pd.DataFrame(np.hstack((X,y)),columns=BDTvars+['y'])
    dfX = df[BDTvars]
    dfY = df['y']

    seed=0
    test_size=0.4
    X_train, X_test, y_train, y_test = train_test_split(dfX, dfY, test_size=test_size, random_state=seed)
    print('Number of events in the training sample:  ', X_train.shape[0])
    print('Number of events in the test sample:  ', X_test.shape[0])
    return BDTvars, df_data_up, df_data_down, df_sig, X_train, X_test, y_train, y_test

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

#PerformGridSearch(X_train, y_train)



def FitModel(X_train, X_test, y_train, y_test):
    model = XGBClassifier(base_score=0.5, booster='gbtree', colsample_bylevel=1,
           colsample_bytree=1.0, gamma=2, learning_rate=0.02,
           max_delta_step=0, max_depth=5, min_child_weight=0.1, missing=None,
           n_estimators=600, n_jobs=1, nthread=1, objective='binary:logistic',
           random_state=0, reg_alpha=0, reg_lambda=1, scale_pos_weight=1,
           seed=None, silent=True, subsample=0.6)

    model.fit(X_train, y_train)

    # Performance on train
    y_pred_train = model.predict(X_train)
    model.fit(X_train, y_train)
    auc_score = roc_auc_score(y_train, y_pred_train)
    print('Using predict function to compute results:')
    print('')
    print("Performance on train sample : ", auc_score)

    # Performance on test
    y_pred_test = model.predict(X_test)
    predictions = [round(value) for value in y_pred_test]
    auc_score = roc_auc_score(y_test, y_pred_test)
    print("Performance on test sample : ", auc_score)

    # evaluate predictions
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy: %.2f%%" % (accuracy * 100.0))
    return model, y_pred_test

def PrintROC(X_test, y_test, y_predicted, model,classifier):
    fpr_lr, tpr_lr, thr = roc_curve(y_test,model.predict_proba(X_test)[:, 1])
    roc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
    areaUC = auc(fpr_lr, tpr_lr)

    print ("Classifier used: ", classifier)
    print (classification_report(y_test, y_predicted,
                                 target_names=["background", "signal"]))
    print ("Area under ROC curve: %.4f"%(roc))
    return fpr_lr, tpr_lr, areaUC

def compare_train_test(model, X_train, y_train, X_test, y_test) :
    """
        Given a trained classifier, get the distribution of the classifier output for the signal and background categ
ories of the training and testing datasets.
        This allows to check for overtraining.
        The output of this function can be passed to other functions for python (python_plot_compare_train_test)
        """

    decisions = []


    for X,y in ((X_train, y_train), (X_test, y_test)) :
        predS = model.predict_proba(X[y>0.5])[:, 1] .ravel()
        predB = model.predict_proba(X[y<0.5])[:, 1] .ravel()

        decisions += [predS, predB]

    # print decisions

    return decisions

def plot_compare_train_test(decisions,bins,classifier):
    """
        Plot the distribution of the classifier output for the signal and background categories of the training and testing datasets.
        This allows to check for overtraining.
        """
    filename='../plots/plt_' + classifier+'_Output.pdf'

    low = min(np.min(d) for d in decisions)
    high = max(np.max(d) for d in decisions)
    #low = 0.
    #high = 1.
    low_high = (low,high)
    # print low_high

    # Plot with python.
    plt.figure()

    plt.hist(decisions[0], color='b', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='S (train)')
    plt.hist(decisions[1], color='r', alpha=0.5, range=low_high, bins=bins, histtype='stepfilled', density=True, label='B (train)')

    hist, bins = np.histogram(decisions[2], bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    width = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.errorbar(center, hist, yerr=err, fmt='o', c='b', label='S (test)')

    hist, bins = np.histogram(decisions[3], bins=bins, range=low_high, density=True)
    scale = len(decisions[2]) / sum(hist)
    err = np.sqrt(hist * scale) / scale

    plt.errorbar(center, hist, yerr=err, fmt='o', c='r', label='B (test)')

    plt.xticks(np.arange(0, 1, step=0.1))
    plt.xlabel("Classifier output")
    plt.ylabel("Arbitrary units")
    plt.legend(loc='best')
    plt.savefig(filename,format='pdf')
    plt.show(block = False)
    return

def PlotSingleROC(fpr, tpr, roc_auc, classifier, color):
    #Plot ROC curve
    plt.figure()
    plt.xlim([-0.01, 1.00])
    plt.ylim([-0.01, 1.01])
    plt.plot(fpr, tpr, lw=1, color = color,
             label='ROC curve (area = %0.3f)' % roc_auc)
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC curve '+ classifier, fontsize=14)
    plt.legend(loc='lower right', fontsize=12)
    plt.savefig('plots/ROCcurve_'+classifier+'.png')
    plt.show(block = False)
    return


def SetBDTvariables():
    #Define BDT variables
    BDTvars =['Lc_FDCHI2_OWNPV','Lc_IPCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT',
              'pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','p_ProbNNp','pi_ProbNNpi','K_ProbNNk']
    return BDTvars

def PlotFeatureImportance(model,classifier):
    features=SetBDTvariables()
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



def PlotLcMass(df_data, BDTcut):
    mLc_s = df_data['Lc_M'].where(df_data['BDT']>BDTcut).dropna()
    mLc_b = df_data['Lc_M'].where(df_data['BDT']<BDTcut).dropna()
    mLc_tot = df_data['Lc_M']

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
    plt_tot.hist(mLc_tot, bins=30, density=False, alpha = 0.5,color = 'green', label='No cut',edgecolor='k')
    plt_tot.hist(mLc_b, bins=30, density=False,alpha = 0.5,color='red', label='BDT < '+str(BDTcut),edgecolor='k')
    plt_tot.hist(mLc_s, bins=30, density=False,alpha = 0.5,color='blue', label='BDT > '+str(BDTcut),edgecolor='k')
    plt_tot.legend(prop={'size': 10})
    plt.savefig('plots/mLcWithBDT.png',bbox_inches='tight')
    plt.show(block = False)


def RunBDT(BDTcut):
    createBDTsamples()
    BDTvars, df_data_up, df_data_down, df_sig, X_train, X_test, y_train, y_test = CreateTrainTest()
    model, y_pred_test=FitModel(X_train, X_test, y_train, y_test)
    print('')
    print('Using predict_proba')
    fpr, tpr, areaUC = PrintROC(X_test, y_test, y_pred_test, model,'xgboost')
    decisions = compare_train_test(model, X_train, y_train, X_test, y_test)
    plot_compare_train_test(decisions,25,'xgboost')
    PlotSingleROC(fpr, tpr, areaUC, 'xgboost', 'navy')
    PlotFeatureImportance(model,'xgboost')
    #To see what are the effects of BDT on Lc mass let's apply the trained model onto the full data dataset
    frames = [df_data_up, df_data_down]
    df_data = pd.concat(frames)
    y_predicted = model.predict_proba(df_data[BDTvars])[:,1]
    df_data = df_data.assign(BDT=y_predicted)
    PlotLcMass(df_data,BDTcut)
    return model




def AddBDTinfo(fname,dtype, model):
    print(fname)
    f = r.TFile(fname,'update')
    if not f:
        print('File NOT found!')
    t = f.Get('DecayTree')
    if not t:
        print('Tree NOT found!')

    bdt = np.zeros(1, dtype=float)
    bdtBranch = t.Branch("bdt",bdt,"bdt/D")

    if dtype!='MC':
        df_data = read_root(fname,'DecayTree', columns=varsData)
    else:
        df_data = read_root(fname,'DecayTree', columns=varsMC)
        df_data =df_data.rename(columns={"p_ProbNNp_corr": "p_ProbNNp","pi_ProbNNpi_corr": "pi_ProbNNpi","K_ProbNNk_corr" : "K_ProbNNk"})

    BDTvars = SetBDTvariables()

    y_predicted = model.predict_proba(df_data[BDTvars])[:, 1]

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        #print(y_predicted[i])
        bdt[0] = y_predicted[i]
        bdtBranch.Fill()

    f.Write()
    f.Close()

