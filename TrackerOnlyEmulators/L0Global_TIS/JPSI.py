import ROOT as r
import math as m

def preselection(t,Lb='Lb'):
    answer = r.kTRUE
    if Lb=='Lb':
            if t.mu_PIDmu<=0: answer = r.kFALSE
            if t.K_PIDK<=0: answer = r.kFALSE
            if t.p_PIDp<=0: answer = r.kFALSE
            if t.Lb_TRUEP_Z<=0: answer = r.kFALSE
            return answer
    if t.mup_PIDmu<=0: answer = r.kFALSE
    if t.mum_PIDmu<=0: answer = r.kFALSE
    if t.mup_PT<=500: answer = r.kFALSE
    if t.mum_PT<=500: answer = r.kFALSE
    if abs(t.Jpsi_MM-3096.9)>=80: answer = r.kFALSE
    if t.Jpsi_ENDVERTEX_CHI2>=16: answer = r.kFALSE

############
# K cuts

    if t.Kplus_IPCHI2_OWNPV<=25: answer = r.kFALSE
    if t.Kplus_PIDK<=0: answer = r.kFALSE
    if t.Kplus_TRACK_GhostProb>=0.35: answer = r.kFALSE

############
# Bp cuts

    if (t.Bplus_ENDVERTEX_CHI2/t.Bplus_ENDVERTEX_NDOF)>=6: answer = r.kFALSE
    if t.Bplus_M<=(5280-50): answer = r.kFALSE
    if t.Bplus_M>=(5280+50): answer = r.kFALSE
    return answer

def TISTOSvsREAL(fjpsi = 'DVntuple_MCMagUp.root',tree = 'TupleBToJpsi_K/DecayTree',var = 'Bplus', tos = 'Bplus_L0MuonDecision_TOS',dataormc = ''):
    cTISTOSvsREAL = r.TCanvas()
    if var=='Bplus' or var=='Bp':
        true =''
        PZ = '_PZ'
        PT = '_PT'
    else:
        true = '_TRUE'
        PZ = 'P_Z'
        PT = 'PT'
    fJPSI_mc = r.TFile(fjpsi)
    tJPSI_mc = fJPSI_mc.Get(tree)
    fJ = r.TFile("L0TIS_efficiency_2D_Julian.root")
    tJPSI_mc.AddFriend("DecayTree",fjpsi[:-5]+"_filtered.root")
    histos = {'ALL':[],'TISTOS':[],'TOS':[],'TIS':[]}
    for iH in histos.keys():
        histos[iH]=fJ.Get("Jpsi_data_eff0").Clone(iH)
        histos[iH].Reset()
        histos[iH].GetXaxis().SetTitleOffset(1)
        histos[iH].GetXaxis().SetTitleOffset(0.5)
        histos[iH].GetYaxis().SetTitleOffset(1)
        histos[iH].GetYaxis().SetTitleOffset(0.5)
        histos[iH].SetDirectory(0)
    print tJPSI_mc.GetEntries()," events in total"
    for i,ie in enumerate(tJPSI_mc):
     if dataormc!='': data = ie.GetLeaf('sweight').GetValue()
     else: data = 1
     if ie.GetLeaf(var+true+PZ).GetValue()>0:
      if preselection(ie,var):
        if i%10000==0: print "processed ",i," events" 
        histos['ALL'].Fill(m.log(ie.GetLeaf(var+true+PZ).GetValue()),m.log(ie.GetLeaf(var+true+PT).GetValue()),data)
        if ie.GetLeaf(tos).GetValue():
            histos['TOS'].Fill(m.log(ie.GetLeaf(var+true+PZ).GetValue()),m.log(ie.GetLeaf(var+true+PT).GetValue()),data)
            if ie.GetLeaf(var+'_L0Global_TIS').GetValue():                                     
                histos['TISTOS'].Fill(m.log(ie.GetLeaf(var+true+PZ).GetValue()),m.log(ie.GetLeaf(var+true+PT).GetValue()),data)
        if ie.GetLeaf(var+'_L0Global_TIS').GetValue():
            histos['TIS'].Fill(m.log(ie.GetLeaf(var+true+PZ).GetValue()),m.log(ie.GetLeaf(var+true+PT).GetValue()),data)
    histos['TIS'].Divide(histos['ALL'])
    histos['TIS'].ProjectionX().GetYaxis().SetRangeUser(2.3,2.8)
    histos['TIS'].ProjectionY().GetYaxis().SetRangeUser(1.,1.8)
    histos['TISTOS'].Divide(histos['TOS'])
    histos['TIS'].SetLineColor(r.kBlue)
    histos['TISTOS'].SetLineColor(r.kGreen)
  #  cTISTOSvsREAL = r.TCanvas()
    cTISTOSvsREAL.Divide(1,3)
    cTISTOSvsREAL.cd(1)
    histos['TIS'].ProjectionX().Draw()
    histos['TISTOS'].ProjectionX().Draw('same')
    cTISTOSvsREAL.cd(2)
    histos['TIS'].ProjectionY().Draw()
    histos['TISTOS'].ProjectionY().Draw('same')
    cTISTOSvsREAL.cd(3)
    hEfficComp = {i:r.TH1F('efficiency_'+i,'efficiency_'+i,32,0,32) for i in ['TIS','TISTOS']}
    for iX in xrange(histos['TIS'].GetXaxis().GetNbins()):
        for iY in xrange(histos['TIS'].GetYaxis().GetNbins()):
            for iT in ['TIS','TISTOS']:
                hEfficComp[iT].SetBinContent(iX*histos['TIS'].GetYaxis().GetNbins()+iY,histos[iT].GetBinContent(iX,iY))
                print iT,"  ",iX*histos['TIS'].GetYaxis().GetNbins()+iY,"   ",histos[iT].GetBinContent(iX,iY)
    hEfficComp['TIS'].SetLineColor(r.kBlue)
    hEfficComp['TISTOS'].SetLineColor(r.kGreen)
    for iH in hEfficComp.keys(): hEfficComp[iH].SetDirectory(0)
    for iH in histos.keys(): histos[iH].SetDirectory(0)
    hEfficComp['TIS'].Draw()
    hEfficComp['TISTOS'].Draw('same')
    cTISTOSvsREAL.g = [histos,hEfficComp]
    cTISTOSvsREAL.SaveAs("TISTOSvsREAL"+fjpsi[-25:-21]+".root")
    cTISTOSvsREAL.SaveAs("TISTOSvsREAL"+fjpsi[-25:-21]+".png")
    histos['TISTOS'].SaveAs("L0TIS_efficiency_2D_2016"+fjpsi[-25:-21]+".root")
    return histos,hEfficComp
    

def makeTISTOSMethod():
    f=r.TFile("DTT_2016_B2JpsiK_RealData_filtered.root")
    t=f.Get("DecayTree")

    mTIS=r.TH2F('TISTOS','TISTOS',160,0,16,200,0,20)
    mTOS=r.TH2F('TOS','TOS',160,0,16,200,0,20) 

    for ie in t:                                                                           
        if ie.Bplus_L0MuonDecision_TOS:                                  
            if ie.Bplus_L0Global_TIS:                                     
                mTIS.Fill(m.log(ie.Bplus_PZ),m.log(ie.Bplus_PT),ie.sweight)


    for ie in t:                                                                           
        if ie.Bplus_L0MuonDecision_TOS:                                  
            mTOS.Fill(m.log(ie.Bplus_PZ),m.log(ie.Bplus_PT),ie.sweight)
        
    co=r.TH2F() 
    mTIS.Copy(co)
    co.Divide(mTOS)
    co.Draw('colz')
