import ROOT as r

f = r.TFile('$FILEDIR/MISID/IntermediateFiles_SS/Templates_Kenr.root','READ')
h_data = f.Get('h_data')
h_muMISID_K = f.Get('h_muMISID_K')
h_muMISID_Pi = f.Get('h_muMISID_Pi')
f1 = r.TFile('$FILEDIR/MISID/IntermediateFiles_SS/muMISIDfractions_Kenr.root','READ')
h_fPi = f1.Get('h_fPi')
h_fK = f1.Get('h_fK')
for ix in range(h_data.GetNbinsX()):
    for iy in range(h_data.GetNbinsY()):
        for iz in range(h_data.GetNbinsZ()):
            print('nbins  = ',ix, iy, iz)
            ndata = h_data.GetBinContent(ix,iy,iz)
            nPi = h_muMISID_Pi.GetBinContent(ix,iy,iz)
            print(ndata, nPi)
            #print(h_data.GetBinContent(ix,iy,iz),h_muMISID_K.GetBinContent(ix,iy,iz),h_muMISID_Pi.GetBinContent(ix,iy,iz), h_fPi.GetBinContent(ix,iy,iz), h_fK.GetBinContent(ix,iy,iz))
