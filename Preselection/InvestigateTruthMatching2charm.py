import ROOT as r

f = r.TFile('/disk/lhcb_data2/RLcMuonic2016/MC_full/Lb_Lc2593Ds_MagUp_full.root')
t = f.Get('tupleout/DecayTree')

t.Draw("mu_MC_MOTHER_ID","Lc_BKGCAT<30 && Lb_BKGCAT<50 && abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13")

t.Draw("mu_MC_MOTHER_ID","Lc_BKGCAT<30 && Lb_BKGCAT<50 && abs(Lb_TRUEID)==5122&&abs(Lc_TRUEID)==4122 &&abs(mu_TRUEID)==13 && (abs(mu_MC_MOTHER_ID)==431 || abs(mu_MC_MOTHER_ID)==421 || abs(mu_MC_MOTHER_ID)==411)")
