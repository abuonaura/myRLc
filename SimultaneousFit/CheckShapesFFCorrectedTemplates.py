import ROOT as r

suffixl = []
MCTO      = False
Iso       = True
Kenr      = False
Lcpi2     = False      

if MCTO:
    MCtype = 'MCTrackerOnly'
else:
    MCtype = 'MCfull'

if Iso: ch = 'Isolated'
if Kenr: ch = 'Kenriched'
if Lcpi2: ch = 'Lcpipi'


def GetFile2Check(ch, MCcat):
    fname = "TemplateFiles/Histos_"+ch+"_"+MCcat+"_LbCorr_FFGstate_FFEstateL_FFEstateH.root"
    print(fname)
    f = r.TFile(fname,'READ')
    return f

def GetFileNoFF(ch,MCcat):
    fname = "TemplateFiles/Histos_"+ch+"_"+MCcat+"_LbCorr_NoFFcorr.root"
    f = r.TFile(fname,'READ')
    return f

def GetMCSampleNames():
    mcsamples = ['Lb_Lcmunu','Lb_Lctaunu','Lb_Lc2625munu','Lb_Lc2625taunu','Lb_Lc2593munu','Lb_Lc2593taunu']
    return mcsamples


def GetHisto(f,sample):
    h = f.Get('h_'+sample)
    return h

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

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

f2Check = GetFile2Check(ch, MCtype)
fNoFF = GetFileNoFF(ch, MCtype)

h2Check = {}
hNoFF = {}

mcsamples = GetMCSampleNames()

'''
for mcsample in mcsamples:
    h2Check[mcsample] = GetHisto(f2Check,mcsample)
    #h2Check[mcsample] = ScaleHisto(h2Check[mcsample],1)
    hNoFF[mcsample] = GetHisto(fNoFF,mcsample)
    #hNoFF[mcsample] = ScaleHisto(hNoFF[mcsample],1)
    c= PlotTemplates(hNoFF[mcsample],'c_'+mcsample)
    c = PlotTemplatesSameCanvas(c,h2Check[mcsample],r.kOrange-3)
    plotdir = 'plots/TemplateFFChecks/FitVars_'+mcsample
    #c.SaveAs(plotdir+'.png')
'''

h = f2Check.Get('h_Lb_Lcmunu')
h1 = fNoFF.Get('h_Lb_Lcmunu')
c1 = PlotTemplates(h1,'c1')
c1 = PlotTemplatesSameCanvas(c1,h,r.kOrange-3)

c2 = PlotTemplates(h,'c2')
c3 = PlotTemplates(h1, 'c3')
