################## - Import
import time
start_time = time.time()
import sys, os
os.environ["CUDA_VISIBLE_DEVICES"] = ""
import zfit
import tensorflow as tf
config = tf.ConfigProto()
config.gpu_options.allow_growth = True
sess = tf.Session(config=config)
zfit.run.sess = sess #run in an instance of RunManager
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams['axes.unicode_minus'] = True
import numpy as np
import uproot
import pandas as pd
zfit.settings.set_seed(100) #set seed
##################

################## - Define function
def Exp_Integral(model, xmin, xmax, obsl):
    def raw_integral(x): return model._unnormalized_pdf(x)/model.params['lambda']
    xpts    = zfit.Data.from_numpy(obs=obsl, array=np.array([xmin, xmax]))
    intgs   = raw_integral(x=xpts)
    return intgs[1] - intgs[0]

def GetData(fname, vars_d, limits):
    t         = uproot.open(fname)["DecayTree"]
    df        = t.pandas.df(vars_d, entrystop=np.inf)
    #Geometrical cuts and Stripping cuts already applied, below trig
    cond   = ((df['Lb_L0Global_TIS'] == 1) | (df['Lc_L0HadronDecision_TOS'] ==1))
    cond   = cond & ((df['Lc_Hlt1TrackMVADecision_TOS'] == 1) | (df['Lc_Hlt1TwoTrackMVADecision_TOS'] ==1))
    cond   = cond & (df['Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS'] ==1)
    #apply muon PID cuts
    cond   = cond & ((df['mu_PIDmu'] > 2) & ((df['mu_PIDmu']-df['mu_PIDK']) > 2) & ((df['mu_PIDmu']-df['mu_PIDp']) > 2) & (df['mu_isMuon']==1))
    #PIDCalib cuts
    cond   = cond & (df['pi_PT'] > 0) & (df['pi_PT'] < 60000) & (df['pi_P'] > 0) & (df['pi_P'] < 200000) 
    cond   = cond & (df['K_PT']  > 0) & (df['K_PT']  < 60000) & (df['K_P']  > 0) & (df['K_P']  < 200000) 
    cond   = cond & (df['p_PT']  > 0) & (df['p_PT']  < 60000) & (df['p_P']  > 0) & (df['p_P']  < 200000) 
    cond   = cond & (df['mu_PT'] > 0) & (df['mu_PT'] < 60000) & (df['mu_P'] > 0) & (df['mu_P'] < 200000) 
    cond   = cond & (df['nTracks'] > 0) & (df['nTracks'] < 700)
    #Offline cuts
    cond   = cond & (df['Lc_M'] > 2230.) & (df['Lc_M'] < 2330.)
    #Below cuts only for FOM optimisation (ie. define bgk region)
    cond   = cond & (df['Lc_M'] > limits[0][0][0]) & (df['Lc_M'] < limits[1][0][0])
    df    = df[cond]
    return df[['Lc_M', 'bdt']]

def GetMC(fname, vars_d, limits):
    t         = uproot.open(fname)["DecayTree"]
    df        = t.pandas.df(vars_d, entrystop=np.inf)
    #Geometrical cuts and Stripping cuts already applied, below trig
    cond   = ((df['Lb_L0Global_TIS'] == 1) | (df['Lc_L0HadronDecision_TOS'] ==1))
    cond   = cond & ((df['Lc_Hlt1TrackMVADecision_TOS'] == 1) | (df['Lc_Hlt1TwoTrackMVADecision_TOS'] ==1))
    cond   = cond & (df['Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS'] ==1)
    #PIDCalib cuts
    cond   = cond & (df['pi_PT'] > 0) & (df['pi_PT'] < 60000) & (df['pi_P'] > 0) & (df['pi_P'] < 200000) 
    cond   = cond & (df['K_PT']  > 0) & (df['K_PT']  < 60000) & (df['K_P']  > 0) & (df['K_P']  < 200000) 
    cond   = cond & (df['p_PT']  > 0) & (df['p_PT']  < 60000) & (df['p_P']  > 0) & (df['p_P']  < 200000) 
    cond   = cond & (df['mu_PT'] > 0) & (df['mu_PT'] < 60000) & (df['mu_P'] > 0) & (df['mu_P'] < 200000) 
    cond   = cond & (df['nTracks'] > 0) & (df['nTracks'] < 700)
    #Offline cuts
    cond   = cond & (df['Lc_M'] > 2230.) & (df['Lc_M'] < 2330.)
    #TruthMatching cuts
    cond   = cond & (abs(df['Lb_TRUEID']) == 5122) & (abs(df['Lc_TRUEID']) == 4122) & (abs(df['p_TRUEID']) == 2212) & (abs(df['K_TRUEID']) == 321) & (abs(df['pi_TRUEID']) == 211) & (abs(df['mu_TRUEID']) == 13)
    cond   = cond & (abs(df['mu_MC_MOTHER_ID']) == 15) & (df['Lc_BKGCAT'] < 30)
    #Below cuts only for FOM optimisation (ie. define signal region)
    cond   = cond & (df['Lc_M'] > limits[0][0][0]) & (df['Lc_M'] < limits[1][0][0])
    num_bef   = np.sum(df['Event_PIDCalibEffWeight'])
    num_after = np.sum(df[cond]['Event_PIDCalibEffWeight'])
    df     = df[cond]
    return df[['Lc_M', 'bdt', 'Event_PIDCalibEffWeight']], num_after, num_bef
###################

##################
#define space to fit in the upperside band region
upper       = (( 2327.,),)
lower       = (( 2310.,),); 
limits      = (lower, upper)
obs         = zfit.Space(obs=['Lc_M'], limits=limits)
#define a space which includes full sig + bkg regs
lowerfull   = (( 2260.,),)
upperfull   = (( 2327.,),)
limitsfull  = (lowerfull, upperfull)
obsfull     = zfit.Space(obs=['Lc_M'], limits=limitsfull)
##################

##################- Get real data read_root
path_data   = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/Data/'
vars_data   = ['Lc_M', 'bdt', 'Lb_L0Global_TIS', 'Lc_L0HadronDecision_TOS', 'Lc_Hlt1TrackMVADecision_TOS', 'Lc_Hlt1TwoTrackMVADecision_TOS', 'Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS', '*_PT', '*_P', 'nTracks', 'mu_PID*', 'mu_isMuon']
dfu = GetData(path_data+"Data_MagUp_2016_MVA.root",   vars_data, limits = limits)
dfd = GetData(path_data+"Data_MagDown_2016_MVA.root", vars_data, limits = limits)
data_df  = pd.concat([dfu,dfd], ignore_index=True)
print('Total Data', data_df.shape)
#plt.hist((zfit.run(data)).ravel(), bins=200)
#plt.savefig('plots/data.pdf')
#plt.clf()
##################

###################
ntot = data_df.shape[0]
params = {}
params['c_combbkg']   = zfit.Parameter(name = 'c_combbkg' ,  value = -0.0015, lower_limit= -1.,      upper_limit=1.,      step_size=0.01, floating=True)
params['yld_combbkg'] = zfit.Parameter(name = 'yld_combbkg', value = ntot,    lower_limit= -2.*ntot, upper_limit=2.*ntot, step_size=5.,   floating=True)
###################

###################
#And make exponential pdf: Defined as exp(param * (x - shift)), the shift is to avoid numerical overflow i.e. exp(large #) gives non-sensical results since 64bit-float cannot handle that mem.
totpdf      = params['yld_combbkg'] * zfit.pdf.Exponential(name = 'Exp_combbkg', obs = obs, lambda_ = params['c_combbkg'])
#define new pdf
totpdffull  = params['yld_combbkg'] * zfit.pdf.Exponential(name = 'Exp_combbkg_full', obs = obsfull, lambda_ = params['c_combbkg'])
totpdffull._set_numerics_data_shift(limits=obsfull)
###################

#################### - set some vals
sppbar_2016     = 144e-6            #https://arxiv.org/pdf/1612.05140.pdf (barns)
f_Lb            = 0.259 * 2 * 0.337 #https://arxiv.org/pdf/1902.06794.pdf, fu~fd (https://arxiv.org/pdf/1306.3663.pdf)
lumi_2016       = 1.67e15           #http://lhcb-public.web.cern.ch/lhcb-public/Images2018/IntRecLumiR12.png (inv barns)
BF_LbToLcmunu   = 6.2e-2; RLc = 0.33
BF_LbToLctaunu  = BF_LbToLcmunu * RLc
BF_TauTomunu    = 17.39e-2
BF_LcTopKpi     = 6.28e-2
####################

#################### - get MC
path_mc           = '/disk/lhcb_data2/amathad/Lb2Lclnu_analysis/MC/'
vars_mc           = ['Lc_M', 'bdt', 'Lb_L0Global_TIS', 'Lc_L0HadronDecision_TOS', 'Lc_Hlt1TrackMVADecision_TOS', 'Lc_Hlt1TwoTrackMVADecision_TOS', 'Lb_Hlt2XcMuXForTauB2XcMuDecision_TOS', 'Lc_BKGCAT', 'Lb_TRUEID', 'Lc_TRUEID', 'p_TRUEID', 'K_TRUEID', 'pi_TRUEID', 'mu_TRUEID', 'mu_MC_MOTHER_ID', '*_PT', '*_P', 'nTracks', 'Event_PIDCalibEffWeight']
dfmcu, nafu, nbfu = GetMC(path_mc+"Lb2Lctaunu_MagDown_2016_PIDCalibAndPIDGen_MVA.root", vars_mc, limits = (lowerfull, lower))
dfmcd, nafd, nbfd = GetMC(path_mc+"Lb2Lctaunu_MagUp_2016_PIDCalibAndPIDGen_MVA.root",   vars_mc, limits = (lowerfull, lower))
mc_df  = pd.concat([dfmcu,dfmcd], ignore_index=True)
print('Total MC', mc_df.shape)
####################

####################
eff_TruthTrigOffl_given_GeomRecoStrip = (nafu+nafd)/(nbfu+nbfd)
print('eff_TruthTrigOffl_given_GeomRecoStrip', eff_TruthTrigOffl_given_GeomRecoStrip)
####################

####################
num_mva_bef                      = np.sum(mc_df['Event_PIDCalibEffWeight'])
N_lbtolcmunu_pt65                = 1e6
num_mva_aft_pt65                 = np.sum(mc_df[mc_df['bdt'] >= 0.65]['Event_PIDCalibEffWeight'])
eff_MVA_given_TruthTrigOffl_pt65 = num_mva_aft_pt65/num_mva_bef #0.014
eff_GeomRecoStrip                = N_lbtolcmunu_pt65/(2. * sppbar_2016 * f_Lb * BF_LbToLcmunu * BF_LcTopKpi * lumi_2016 * eff_TruthTrigOffl_given_GeomRecoStrip * eff_MVA_given_TruthTrigOffl_pt65)
print('eff_GeomRecoStrip', eff_GeomRecoStrip)
####################

mvacuts = np.linspace(0.0,0.9,300)[:-1]
for mvacutval in mvacuts:
    print('MVA cut', mvacutval)

    ############
    #make data with mvacuts and set yld value
    data     = zfit.Data.from_pandas(data_df[data_df['bdt'] > mvacutval]['Lc_M'], obs=obs, weights=None)
    ntot     = float(((zfit.run(data)).ravel()).shape[0])
    params['yld_combbkg'].set_value(ntot)
    print('Before Fit', zfit.run(params['yld_combbkg']), zfit.run(params['c_combbkg']))
    ############

    ############
    #declare the nll and the minimiser and conduct the fit
    nll         = zfit.loss.ExtendedUnbinnedNLL(model=totpdf, data=data, fit_range=obs, constraints=None)
    minimizer   = zfit.minimize.Minuit() #verbosity=10
    cond = False
    nfits = 0
    while (not cond) and nfits < 10:
        print('Conducting nfits', nfits)
        for param in list(params.values()): print('Before Randomisation', param.name, zfit.run(param), zfit.run(param.lower_limit), zfit.run(param.upper_limit))
        for param in list(params.values()): 
            if 'yld' not in param.name: 
                param.randomize(minval=zfit.run(param.lower_limit), maxval=zfit.run(param.upper_limit)) #randomize
    
        for param in list(params.values()): print('After Randomisation', param.name, zfit.run(param))
        result = minimizer.minimize(loss=nll)
        rinfo  = result.info['original']
        #print(result.info)
        #cond = (rinfo['is_valid'] and rinfo['has_valid_parameters'] and rinfo['has_covariance'] and not rinfo['hesse_failed'] and not rinfo['is_above_max_edm'] and not rinfo['has_reached_call_limit'])
        cond = (rinfo['is_valid'] and rinfo['has_valid_parameters'] and rinfo['has_covariance'] and not rinfo['hesse_failed'])
        print("Has a good fit been found" , cond)
        if cond:
            for p, errors in result.hesse().items(): 
                print(p.name, result.params[p]['value'], '+/-', errors['error'])
                params[p.name].set_value(result.params[p]['value'])
        nfits += 1
    #tf.get_default_graph().finalize() #need to do this to use the same default graph
    ##################

    #################### - plot the fit result
    #his1, bins1, _  = plt.hist((zfit.run(data)).ravel(), bins=200)
    #pts             = bins1[:-1] + (bins1[1:] - bins1[:-1])/2.
    #p_bkg           = zfit.run(params['yld_combbkg']) * zfit.run(totpdf.pdf(pts))/(zfit.run(totpdf.pdf(pts))).sum()
    #plt.plot(pts, p_bkg)
    #plt.savefig('plots/fit.pdf')
    #####################

    ############
    #calculate num of bkg cands in signal reg.
    print('After Fit', zfit.run(params['yld_combbkg']), zfit.run(params['c_combbkg']))
    intg_bkg_reg = zfit.run(Exp_Integral(model = totpdffull, xmin = lower[0][0],     xmax = upper[0][0],  obsl = obsfull))
    intg_sig_reg = zfit.run(Exp_Integral(model = totpdffull, xmin = lowerfull[0][0], xmax = lower[0][0],  obsl = obsfull))
    B            = ntot * intg_sig_reg/intg_bkg_reg
    ############

    ####################
    num_mva_aft = np.sum(mc_df[mc_df['bdt'] >= mvacutval]['Event_PIDCalibEffWeight'])
    eff_MVA_given_TruthTrigOffl     = num_mva_aft/num_mva_bef
    print('eff_MVA_given_TruthTrigOffl', eff_MVA_given_TruthTrigOffl)
    ####################

    ####################
    toteff          = eff_GeomRecoStrip * eff_TruthTrigOffl_given_GeomRecoStrip * eff_MVA_given_TruthTrigOffl
    S               = 2. * sppbar_2016 * f_Lb * lumi_2016 * BF_LbToLctaunu * BF_LcTopKpi * BF_TauTomunu * toteff
    ####################

    ####################
    FoM1 = S/np.sqrt(S+B)
    FoM2 = S/(S+B) * FoM1
    ####################

    ####################
    print(mvacutval,',',S,',',B,',',eff_MVA_given_TruthTrigOffl,',',FoM1,',',FoM2)
    if not os.path.isdir('FOMcourse'): os.mkdir('FOMcourse')
    f=open("FOMcourse/FOM.csv", "a+")
    f.write(str(mvacutval)+','+str(S)+','+str(B)+','+str(eff_MVA_given_TruthTrigOffl)+','+str(FoM1)+','+str(FoM2)+'\n')
    f.close()
    ####################
