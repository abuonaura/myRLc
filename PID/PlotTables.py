import ROOT as r

def GetTable(fname, histname):
    f = r.TFile(fname,'READ')
    h = f.Get(histname)
    h.SetDirectory(0)
    return h

def GetTableProjections(h):
    return [h.ProjectionX(), h.ProjectionY(), h.ProjectionZ()]

def PlotTable(name, h):
    c = r.TCanvas('c_'+name,'c_'+name, 1500,500)
    c.Divide(3,1)
    c.cd(1)
    h[0].Draw('hist')
    c.cd(2)
    h[1].Draw('hist')
    c.cd(3)
    h[2].Draw('hist')
    return c

particles = ['K','Pi']
cuts = ['IsMuon==0 && DLLK>4.0', 'IsMuon==0 && DLLK>4.0 && DLLp-DLLK<0', 'IsMuon==0 && DLLK>4.0 && DLLp-DLLK>=0']
folder = '/disk/lhcb_data2/RLcMuonic2016/HistoPID/KenrichedSelection/'
polarities = ['MagUp','MagDown']

h,c = {},{}
for particle in particles:
    h[particle],c[particle] = {},{}
    for polarity in polarities:
        h[particle][polarity],c[particle][polarity] = {},{}
        fname = folder+'PerfHists_'+particle+'_Turbo16_'+polarity+'_Brunel_P_Brunel_ETA_nTracks_Brunel.root'
        for cut in cuts:
            hname = particle+'_'+cut+'_All'
            h[particle][polarity][cut] = GetTable(fname, hname)
            print(h[particle][polarity][cut])
            hpj = GetTableProjections(h[particle][polarity][cut])
            c[particle][polarity][cut] = PlotTable(particle+'_'+polarity+'_'+cut,hpj)
            c[particle][polarity][cut].Draw()
