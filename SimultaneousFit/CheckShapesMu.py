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
    shades = {'blue':600,'orange':797,'red':634,'magenta':617,'black':1}
    h1.SetLineColor(shades[colour])
    h1.SetMarkerColor(shades[colour])
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

h_mu = GetTemplate('mu')
h_FFmu = GetFFTemplate('mu_mean')
h_FFmu_max = GetFFTemplate('mu_max')
h_FFmu_min = GetFFTemplate('mu_min')

r.gStyle.SetOptStat(0)

c = r.TCanvas('c','c: Comparison mu before/after FF correction',1500,500)
h_mu = AssignColour(h_mu,'blue')
h_FFmu = AssignColour(h_FFmu, 'orange')
h_mu_q2, h_mu_El, h_mu_Mm= GetProjections(h_mu)
h_FFmu_q2, h_FFmu_El, h_FFmu_Mm= GetProjections(h_FFmu)
c.Divide(3,1)
c.cd(1)
h_mu_q2.Draw('hist')
h_FFmu_q2.Draw('hist same')
c.cd(2)
h_mu_El.Draw('hist')
h_FFmu_El.Draw('hist same')
c.cd(3)
h_mu_Mm.Draw('hist')
h_FFmu_Mm.Draw('hist same')
l = r.TLegend(0.1,0.7,0.4,0.9)
l.AddEntry(h_mu_Mm,"Before FF corr",'l')
l.AddEntry(h_FFmu_Mm,"After FF corr",'l')
l.SetTextSize(0.04)
l.Draw()

c1 = r.TCanvas('c1','c1: Comparison shapes with FF correction',1500,500)
h_FFmu_min = AssignColour(h_FFmu_min,'blue')
h_FFmu = AssignColour(h_FFmu, 'black')
h_FFmu_max = AssignColour(h_FFmu_max, 'red')
h_FFmu_min_q2, h_FFmu_min_El, h_FFmu_min_Mm = GetProjections(h_FFmu_min)
h_FFmu_q2, h_FFmu_El, h_FFmu_Mm = GetProjections(h_FFmu)
h_FFmu_max_q2, h_FFmu_max_El, h_FFmu_max_Mm = GetProjections(h_FFmu_max)
c1.Divide(3,1)
c1.cd(1)
h_FFmu_max_q2.Draw('PE0')
h_FFmu_q2.Draw('PE0 same')
h_FFmu_min_q2.Draw('PE0 same')
c1.cd(2)
h_FFmu_max_El.Draw('PE0')
h_FFmu_El.Draw('PE0 same')
h_FFmu_min_El.Draw('PE0 same')
c1.cd(3)
h_FFmu_max_Mm.Draw('PE0')
h_FFmu_Mm.Draw('PE0 same')
h_FFmu_min_Mm.Draw('PE0 same')
l1 = r.TLegend(0.1,0.7,0.4,0.9)
l1.AddEntry(h_FFmu_min_Mm,"Min",'lp')
l1.AddEntry(h_FFmu_Mm,"Nominal",'lp')
l1.AddEntry(h_FFmu_max_Mm,'Max','lp')
l1.SetTextSize(0.04)
l1.Draw()
