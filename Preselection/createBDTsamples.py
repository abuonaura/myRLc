import ROOT as r


datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/'
files = {'bkg':{'MagUp':datadir+'Data/Lb_Data_MagUp_reduced.root','MagDown':datadir+'Data/Lb_Data_MagDown_reduced.root'},'sig':datadir+'MC/Lb_Lctaunu_PID_reduced.root'}

varsMC = ['Lc_M','Lc_BKGCAT','Lc_FDCHI2_OWNPV','Lc_IPCHI2_OWNPV','Lc_ENDVERTEX_CHI2','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','p_ProbNNp_corr','pi_ProbNNpi_corr','K_ProbNNk_corr']
f_sig = r.TFile(files['sig'],'READ')
t_sig = f_sig.Get('DecayTree')
t_sig.SetBranchStatus('*',0)
for var in varsMC:
    t_sig.SetBranchStatus(var,1)
cut_sig ='Lc_M>2260. && Lc_M<2310. && Lc_BKGCAT<30.'


ofname_sig =files['sig'][0:-5]+'_forBDT.root'
of_sig = r.TFile(ofname_sig,'recreate')
ot_sig = r.TTree('DecayTree','DecayTree')
print('>>> Copying signal file: '+files['sig'])
ot_sig =t_sig.CopyTree(cut_sig)
print('>>> Created signal file for BDT: '+ofname_sig)
ot_sig.Write()
of_sig.Close()
of_sig.Close()


varsData = ['Lc_M','Lc_FDCHI2_OWNPV','Lc_ENDVERTEX_CHI2','Lc_IPCHI2_OWNPV','pi_PT','p_PT','K_PT','pi_MINIPCHI2','p_MINIPCHI2','K_MINIPCHI2','p_ProbNNp','pi_ProbNNpi','K_ProbNNk']
f_bkg_MagUp = r.TFile(files['bkg']['MagUp'],'READ')
t_bkg_MagUp = f_bkg_MagUp.Get('DecayTree')
t_bkg_MagUp.SetBranchStatus('*',0)
for var in varsData:
    t_bkg_MagUp.SetBranchStatus(var,1)
cut_bkg = 'Lc_M<2260. || Lc_M>2310.'

ofname_bkg_MagUp =files['bkg']['MagUp'][0:-5]+'_forBDT.root'
ofname_bkg_MagDown =files['bkg']['MagDown'][0:-5]+'_forBDT.root'
of_bkg_MagUp = r.TFile(ofname_bkg_MagUp,'recreate')
ot_bkg_MagUp = r.TTree('DecayTree','DecayTree')
print('>>> Copying bkg file: '+files['bkg']['MagUp'])
ot_bkg_MagUp =t_bkg_MagUp.CopyTree(cut_bkg,'',100000)
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
ot_bkg_MagDown =t_bkg_MagDown.CopyTree(cut_bkg,'',100000)
print('>>> Created bkg file for BDT: '+ofname_bkg_MagDown)
ot_bkg_MagDown.Write()
of_bkg_MagDown.Close()
f_bkg_MagDown.Close()


