import ROOT as r

datadir = '/disk/gangadir/buonaura/gangadir/workspace/buonaura/LocalXML/Datasets/MC/'

f = r.TFile(datadir+'Lb_LcDs_PID_reduced_preselected.root','READ')
t = f.Get('DecayTree')

cut2body = 'Lb_TrueHadron_D0_ID!=0&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID==0'
cut3body = 'Lb_TrueHadron_D0_ID!=0&&Lb_TrueHadron_D1_ID!=0&&Lb_TrueHadron_D2_ID!=0'
TruthMatching = 'Lb_BKGCAT<50&&Lc_BKGCAT<30'
cutIsoBDT = 'Lb_ISOLATION_BDT>0.35'
fullIsolation = '(Lb_ISOLATION_BDT>0.35&&Lb_ISOLATION_PIDK>4.)||(Lb_ISOLATION_BDT2>0.35&&Lb_ISOLATION_PIDK2>4.)'


h_2bd = r.TH1F('h_2bd','h_2bd',50,-2,1)
h_3bd = r.TH1F('h_3bd','h_3bd',50,-2,1)
h_2bd_isocut = r.TH1F('h_2bd_isocut','h_2bd_isocut',50,-2,1)
h_3bd_isocut = r.TH1F('h_3bd_isocut','h_3bd_isocut',50,-2,1)
h_2bd_fulliso = r.TH1F('h_2bd_fulliso','h_2bd_isocut',50,-2,1)
h_3bd_fulliso = r.TH1F('h_3bd_fulliso','h_3bd_isocut',50,-2,1)


t.Draw('Lb_ISOLATION_BDT>>h_2bd',cut2body+'&&'+TruthMatching)
h_2bd = r.gPad.GetPrimitive("h_2bd")
t.Draw('Lb_ISOLATION_BDT>>h_2bd_isocut',cut2body+'&&'+TruthMatching+'&&'+cutIsoBDT)
h_2bd_isocut = r.gPad.GetPrimitive("h_2bd_isocut")
t.Draw('Lb_ISOLATION_BDT>>h_2bd_fulliso',cut2body+'&&'+TruthMatching+'&&'+fullIsolation)
h_2bd_fulliso = r.gPad.GetPrimitive("h_2bd_fulliso")

t.Draw('Lb_ISOLATION_BDT>>h_3bd',cut3body+'&&'+TruthMatching)
h_3bd = r.gPad.GetPrimitive("h_3bd")
t.Draw('Lb_ISOLATION_BDT>>h_3bd_isocut',cut3body+'&&'+TruthMatching+'&&'+cutIsoBDT)
h_3bd_isocut = r.gPad.GetPrimitive("h_3bd_isocut")
t.Draw('Lb_ISOLATION_BDT>>h_3bd_fulliso',cut3body+'&&'+TruthMatching+'&&'+fullIsolation)
h_3bd_fulliso = r.gPad.GetPrimitive("h_3bd_fulliso")

c = r.TCanvas('c','2body',500,500)
h_2bd.Draw()
h_2bd_isocut.SetFillColor(r.kBlue)
h_2bd_isocut.Draw('same')
#h_2bd_fulliso.SetFillStyle(3004)
#h_2bd_fulliso.SetFillColor(r.kBlue)
#h_2bd_fulliso.Draw('same')

c1 = r.TCanvas('c1','3body',500,500)
h_3bd.Draw()
h_3bd_isocut.SetFillColor(r.kRed)
h_3bd_isocut.Draw('same')

print('fraction of 2body and 3 body decays selected by ', cutIsoBDT)
print(' - 2 body: {:.4f}'.format(h_2bd_isocut.GetEntries()*1./h_2bd.GetEntries()))
print(' - 3 body: {:.4f}'.format(h_3bd_isocut.GtgetEntries()*1./h_3bd.GetEntries()))
print('')
print('fraction of 2body and 3 body decays selected by ', fullIsolation)
print(' - 2 body: {:.4f}'.format(h_2bd_fulliso.GetEntries()*1./h_2bd.GetEntries()))
print(' - 3 body: {:.4f}'.format(h_3bd_fulliso.GetEntries()*1./h_3bd.GetEntries()))


#h_3bd_fulliso.SetFillStyle(3004)
#h_3bd_fulliso.SetFillColor(r.kRed)
#h_3bd_fulliso.Draw('same')


