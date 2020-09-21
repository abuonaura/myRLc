import ROOT as r

def GetTemplate(name):
    f = r.TFile('RootFiles/Histos_Isolated_MCfull.root','READ')
    h = f.Get('h_'+name)
    h.SetDirectory(0)
    return h

def GetFFTemplate(name):
    f = r.TFile('RootFiles/DemoHistosLattice_Isolated_MCfull.root','READ')
    h = f.Get('h_w_'+name)
    h.SetDirectory(0)
    return h

def AssignColour(h, colour):
    h1 = h.Clone(h.GetName())
    h1.SetDirectory(0)
    shades = {'blue':600,'orange':797,'red':634,'magenta':617,'black':1,'violet':877}
    h1.SetLineColor(shades[colour])
    h1.SetMarkerColor(shades[colour])
    h1.SetMarkerStyle(8)
    return h1

def GetProjections(h):
    hx = h.ProjectionX()
    hx.GetXaxis().SetTitle('q^{2}')
    hx.SetDirectory(0)
    hy = h.ProjectionY()
    hy.GetXaxis().SetTitle('E_{l}')
    hy.SetDirectory(0)
    hz = h.ProjectionZ()
    hz.GetXaxis().SetTitle('M^{2}_{miss}')
    hz.SetDirectory(0)
    return hx, hy, hz

def ScaleHisto(h,value):
    h1 = h.Clone(h.GetName()+'_norm')
    scale = value*1./h1.Integral()
    h1.Scale(scale)
    return h1

def GetProjectionsNormalised(h, value):
    hx, hy, hz = GetProjections(h)
    hx_norm = ScaleHisto(hx,value)
    hx_norm.SetDirectory(0)
    hy_norm = ScaleHisto(hy,value)
    hy_norm.SetDirectory(0)
    hz_norm = ScaleHisto(hz,value)
    hz_norm.SetDirectory(0)
    return hx_norm, hy_norm, hz_norm


h_tau = GetTemplate('tau')
h_FFtau = GetFFTemplate('tau_mean')
h_FFtau_max = GetFFTemplate('tau_max')
h_FFtau_min = GetFFTemplate('tau_min')

r.gStyle.SetOptStat(0)

c = r.TCanvas('c','c: Comparison tau before/after FF correction',1500,500)
h_tau = AssignColour(h_tau,'blue')
h_FFtau = AssignColour(h_FFtau, 'orange')
h_tau_q2, h_tau_El, h_tau_Mm= GetProjections(h_tau)
h_FFtau_q2, h_FFtau_El, h_FFtau_Mm= GetProjections(h_FFtau)
c.Divide(3,1)
c.cd(1)
h_tau_q2.Draw('hist')
h_FFtau_q2.Draw('hist same')
c.cd(2)
h_tau_El.Draw('hist')
h_FFtau_El.Draw('hist same')
c.cd(3)
h_tau_Mm.Draw('hist')
h_FFtau_Mm.Draw('hist same')
l = r.TLegend(0.1,0.7,0.4,0.9)
l.AddEntry(h_tau_Mm,"Before FF corr",'l')
l.AddEntry(h_FFtau_Mm,"After FF corr",'l')
l.SetTextSize(0.04)
l.Draw()

c1 = r.TCanvas('c1','c1: Comparison shapes with FF correction',1500,500)
h_FFtau_min = AssignColour(h_FFtau_min,'blue')
h_FFtau = AssignColour(h_FFtau, 'black')
h_FFtau_max = AssignColour(h_FFtau_max, 'red')
h_FFtau_min_q2, h_FFtau_min_El, h_FFtau_min_Mm = GetProjections(h_FFtau_min)
h_FFtau_q2, h_FFtau_El, h_FFtau_Mm = GetProjections(h_FFtau)
h_FFtau_max_q2, h_FFtau_max_El, h_FFtau_max_Mm = GetProjections(h_FFtau_max)
c1.Divide(3,1)
c1.cd(1)
h_FFtau_max_q2.Draw('PE0')
h_FFtau_q2.Draw('PE0 same')
h_FFtau_min_q2.Draw('PE0 same')
c1.cd(2)
h_FFtau_max_El.Draw('PE0')
h_FFtau_El.Draw('PE0 same')
h_FFtau_min_El.Draw('PE0 same')
c1.cd(3)
h_FFtau_max_Mm.Draw('PE0')
h_FFtau_Mm.Draw('PE0 same')
h_FFtau_min_Mm.Draw('PE0 same')
l1 = r.TLegend(0.1,0.7,0.4,0.9)
l1.AddEntry(h_FFtau_min_Mm,"Min",'lp')
l1.AddEntry(h_FFtau_Mm,"Nominal",'lp')
l1.AddEntry(h_FFtau_max_Mm,'Max','lp')
l1.SetTextSize(0.04)
l1.Draw()

h_MISID = GetTemplate('MISID')
h_Comb = GetTemplate('Combinatorial')
h_mu = GetTemplate('mu')
h_startau = GetTemplate('startau')
h_starmu = GetTemplate('starmu')

c2 = r.TCanvas('c2','c2: Comparison tau/MISID',1500,500)
h_MISID = AssignColour(h_MISID, 'orange')
h_tau_q2_norm, h_tau_El_norm, h_tau_Mm_norm = GetProjectionsNormalised(h_tau,1.)
h_MISID_q2_norm, h_MISID_El_norm, h_MISID_Mm_norm = GetProjectionsNormalised(h_MISID,1.)
c2.Divide(3,1)
c2.cd(1)
h_tau_q2_norm.Draw('hist')
h_MISID_q2_norm.Draw('hist same')
c2.cd(2)
h_tau_El_norm.Draw('hist')
h_MISID_El_norm.Draw('hist same')
c2.cd(3)
h_tau_Mm_norm.Draw('hist')
h_MISID_Mm_norm.Draw('hist same')
l2 = r.TLegend(0.7,0.7,0.9,0.9)
l2.AddEntry(h_tau_Mm_norm,"#tau",'l')
l2.AddEntry(h_MISID_Mm_norm,"MISID",'l')
l2.SetTextSize(0.04)
l2.Draw()

c3 = r.TCanvas('c3','c3: Comparison tau/Comb',1500,500)
h_Comb = AssignColour(h_Comb, 'orange')
h_Comb_q2_norm, h_Comb_El_norm, h_Comb_Mm_norm = GetProjectionsNormalised(h_Comb,1.)
c3.Divide(3,1)
c3.cd(1)
h_tau_q2_norm.Draw('hist')
h_Comb_q2_norm.Draw('hist same')
c3.cd(2)
h_tau_El_norm.Draw('hist')
h_Comb_El_norm.Draw('hist same')
c3.cd(3)
h_tau_Mm_norm.Draw('hist')
h_Comb_Mm_norm.Draw('hist same')
l3 = r.TLegend(0.7,0.7,0.9,0.9)
l3.AddEntry(h_tau_Mm_norm,"#tau",'l')
l3.AddEntry(h_Comb_Mm_norm,"Comb",'l')
l3.SetTextSize(0.04)
l3.Draw()

c4 = r.TCanvas('c4','c4: Comparison tau/mu',1500,500)
h_mu = AssignColour(h_mu, 'orange')
h_mu_q2_norm, h_mu_El_norm, h_mu_Mm_norm = GetProjectionsNormalised(h_mu,1.)
c4.Divide(3,1)
c4.cd(1)
h_tau_q2_norm.Draw('hist')
h_mu_q2_norm.Draw('hist same')
c4.cd(2)
h_tau_El_norm.Draw('hist')
h_mu_El_norm.Draw('hist same')
c4.cd(3)
h_tau_Mm_norm.Draw('hist')
h_mu_Mm_norm.Draw('hist same')
l4 = r.TLegend(0.7,0.7,0.9,0.9)
l4.AddEntry(h_tau_Mm_norm,"#tau",'l')
l4.AddEntry(h_mu_Mm_norm,"mu",'l')
l4.SetTextSize(0.04)
l4.Draw()

c5 = r.TCanvas('c5','c5: Comparison tau/startau',1500,500)
h_startau = AssignColour(h_startau, 'orange')
h_startau_q2_norm, h_startau_El_norm, h_startau_Mm_norm = GetProjectionsNormalised(h_startau,1.)
c5.Divide(3,1)
c5.cd(1)
h_tau_q2_norm.Draw('hist')
h_startau_q2_norm.Draw('hist same')
c5.cd(2)
h_tau_El_norm.Draw('hist')
h_startau_El_norm.Draw('hist same')
c5.cd(3)
h_tau_Mm_norm.Draw('hist')
h_startau_Mm_norm.Draw('hist same')
l5 = r.TLegend(0.7,0.7,0.9,0.9)
l5.AddEntry(h_tau_Mm_norm,"#tau",'l')
l5.AddEntry(h_startau_Mm_norm,"startau",'l')
l5.SetTextSize(0.04)
l5.Draw()

c6 = r.TCanvas('c6','c6: Comparison tau/starmu',1500,500)
h_starmu = AssignColour(h_starmu, 'orange')
h_starmu_q2_norm, h_starmu_El_norm, h_starmu_Mm_norm = GetProjectionsNormalised(h_starmu,1.)
c6.Divide(3,1)
c6.cd(1)
h_tau_q2_norm.Draw('hist')
h_starmu_q2_norm.Draw('hist same')
c6.cd(2)
h_tau_El_norm.Draw('hist')
h_starmu_El_norm.Draw('hist same')
c6.cd(3)
h_tau_Mm_norm.Draw('hist')
h_starmu_Mm_norm.Draw('hist same')
l6 = r.TLegend(0.7,0.7,0.9,0.9)
l6.AddEntry(h_tau_Mm_norm,"#tau",'l')
l6.AddEntry(h_starmu_Mm_norm,"starmu",'l')
l6.SetTextSize(0.04)
l6.Draw()

h_MISID = GetTemplate('MISID')
h_Comb = GetTemplate('Combinatorial')
h_mu = GetTemplate('mu')
h_startau = GetTemplate('startau')
h_starmu = GetTemplate('starmu')

c7 = r.TCanvas('c7','c7: all together', 1500,500)
c7.Divide(3,1)
h_mu = AssignColour(h_mu,'orange')
h_startau = AssignColour(h_startau,'black')
h_starmu = AssignColour(h_starmu,'red')
h_MISID = AssignColour(h_MISID,'magenta')
h_Comb = AssignColour(h_Comb,'violet')
h_tau1_q2_norm, h_tau1_El_norm, h_tau1_Mm_norm = GetProjectionsNormalised(h_tau,1.)
h_mu1_q2_norm, h_mu1_El_norm, h_mu1_Mm_norm = GetProjectionsNormalised(h_mu,1.)
h_startau1_q2_norm, h_startau1_El_norm, h_startau1_Mm_norm = GetProjectionsNormalised(h_startau,1.)
h_starmu1_q2_norm, h_starmu1_El_norm, h_starmu1_Mm_norm = GetProjectionsNormalised(h_starmu,1.)
h_MISID1_q2_norm, h_MISID1_El_norm, h_MISID1_Mm_norm = GetProjectionsNormalised(h_MISID,1.)
h_Comb1_q2_norm, h_Comb1_El_norm, h_Comb1_Mm_norm = GetProjectionsNormalised(h_Comb,1.)
c7.cd(1)
h_starmu1_q2_norm.GetYaxis().SetRangeUser(0,0.8)
h_starmu1_q2_norm.Draw('hist')
h_mu1_q2_norm.Draw('hist same')
h_tau1_q2_norm.Draw('PE0 same')
h_startau1_q2_norm.Draw('PE0 same')
h_MISID1_q2_norm.Draw('hist same')
h_Comb1_q2_norm.Draw('hist same')
c7.cd(2)
h_starmu1_El_norm.GetYaxis().SetRangeUser(0,0.3)
h_starmu1_El_norm.Draw('hist')
h_mu1_El_norm.Draw('hist same')
h_tau1_El_norm.Draw('P0 same')
h_startau1_El_norm.Draw('PE0 same')
h_MISID1_El_norm.Draw('hist same')
h_Comb1_El_norm.Draw('hist same')
c7.cd(3)
h_starmu1_Mm_norm.GetYaxis().SetRangeUser(0,0.7)
h_starmu1_Mm_norm.Draw('hist')
h_mu1_Mm_norm.Draw('hist same')
h_tau1_Mm_norm.Draw('PE0 same')
h_startau1_Mm_norm.Draw('PE0 same')
h_MISID1_Mm_norm.Draw('hist same')
h_Comb1_Mm_norm.Draw('hist same')
l7 = r.TLegend(0.7,0.7,0.9,0.9)
l7.AddEntry(h_tau1_Mm_norm,'#tau','lp')
l7.AddEntry(h_mu1_Mm_norm,'mu','l')
l7.AddEntry(h_startau1_Mm_norm,'startau','lp')
l7.AddEntry(h_starmu1_Mm_norm,'starmu','l')
l7.AddEntry(h_MISID1_Mm_norm,'MISID','l')
l7.AddEntry(h_Comb1_Mm_norm,'Comb','l')
l7.Draw()