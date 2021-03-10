import ROOT as r


def PlotTemplates(h,name):
    h.SetLineColor(r.kAzure-3)
    c = r.TCanvas(name,'',1500,500)
    c.Divide(3,1)
    c.cd(1)
    h.ProjectionX().Draw('hist')
    c.cd(2)
    h.ProjectionY().Draw('hist')
    c.cd(3)
    h.ProjectionZ().Draw('hist')
    return c

def PlotTemplatesSameCanvas(c,h,color):
    h.SetLineColor(color)
    c.cd(1)
    h.ProjectionX().Draw('hist same')
    c.cd(2)
    h.ProjectionY().Draw('hist same')
    c.cd(3)
    h.ProjectionZ().Draw('hist same')
    return c

categories = ['Isolated','Kenriched','Lcpipi']

for category in categories:
    fnew = r.TFile('/home/hep/buonaura/Analysis/RLc/scripts/SimultaneousFit/RootFiles/Histos_'+category+'_MCfull.root')
    hmu = fnew.Get('h_mu')
    htau = fnew.Get('h_tau')
    fold = r.TFile('/home/hep/buonaura/Analysis/RLc/scripts/SimultaneousFit/RootFiles/DemoHistosLattice_'+category+'_MCfull.root')
    hmuOLD = fold.Get('h_w_mu_mean')
    htauOLD = fold.Get('h_w_tau_mean')
    
    c = PlotTemplates(hmu,'c')
    c = PlotTemplatesSameCanvas(c,hmuOLD,r.kOrange-3)
    c.SaveAs('plots/ComparisonFFcorrections_mu_'+category+'.png')
    c1 = PlotTemplates(htau,'c1')
    c1 = PlotTemplatesSameCanvas(c1,htauOLD,r.kOrange-3)
    c1.SaveAs('plots/ComparisonFFcorrections_tau_'+category+'.png')
