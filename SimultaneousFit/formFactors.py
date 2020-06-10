import ROOT as r
from math import *
import multi

m_B = 5.279
m_pi = 0.1348
m_Bc = 6.276
m_Lc = 2.28646
m_Lb = 5.6195
Gf = 1.1663787E-5
m_m = 0.106
m_t = 1.77
pi = 3.14
ei = 0
V = 999999

def obtaineTheorieParameters():
    multiGauss = multi.PyMultiGauss()
    diffParam = multiGauss.get_parameters()
    return diffParam

# function is corresponfing to the equition (72)
def defZ(q2,tp,t0):
	z = (sqrt(tp-q2)-sqrt(tp-t0))/(sqrt(tp-q2)+sqrt(tp-t0))
	return z

def defParamGener(generPoints, typ,count=0):
        if typ=='fp' or typ=='ft': deltaF = 0.056
        if typ=='f0': deltaF = 0.449
        if typ=='gp' or typ=='gt': deltaF = 0.492
        if typ=='g0': deltaF = 0
        for iK in generPoints.keys():
                    if iK=='xa0'+typ: a0 = generPoints[iK][count]
                    if iK=='xa1'+typ: a1 = generPoints[iK][count]
        if typ=="gt": a0 = generPoints["xa0gp"][count]
        m_pole = m_Bc + deltaF
        tp = m_pole*m_pole
        t0 = (m_Lb - m_Lc)*(m_Lb - m_Lc)        # (73)
        qmax2 = (m_Lb - m_Lc)*(m_Lb - m_Lc)
        return [m_pole,a0,a1,tp,t0,qmax2]
	
def defParam(typ):
	if typ=='fp' or typ=='ft': deltaF = 0.056
	if typ=='f0': deltaF = 0.449
	if typ=='gp' or typ=='gt': deltaF = 0.492
	if typ=='g0': deltaF = 0
	if typ=='fp':
		a0 = 0.8146
		a1 = -4.899
        if typ=='ft':
                a0 = 1.0780
                a1 = -6.4170
        if typ=='f0':
                a0 = 0.7439
                a1 = - 4.6480
        if typ=='gp':
                a0 = 0.6847
                a1 = -4.4310
        if typ=='gt':
                a0 = 0.6847
                a1 = -4.4630
        if typ=='g0':
                a0 = 0.7396
                a1 = -4.3660
	m_pole = m_Bc + deltaF
	tp = m_pole*m_pole
	t0 = (m_Lb - m_Lc)*(m_Lb - m_Lc) 	# (73)
	qmax2 = (m_Lb - m_Lc)*(m_Lb - m_Lc)
	return [m_pole,a0,a1,tp,t0,qmax2]

# function is corresponfing to the equition (79)
def count_f(q2,m_pole,a0,a1,z):
	return (1/(1-q2/(m_pole*m_pole)))*(a0+a1*z)
def all_f(q2,generPoints,count=0,multiGauss=1):
	formFack = {}
	fptypes=['fp','ft','f0','gp','gt','g0']
	for iT in fptypes:
                if multiGauss:m_pole,a0,a1,tp,t0,qmax2 = defParamGener(generPoints, iT,count)
		else: m_pole,a0,a1,tp,t0,qmax2 = defParam(iT)
		z = defZ(q2,tp,t0)
		formFack[iT] = count_f(q2,m_pole,a0,a1,z)
	return formFack

def partWidth(q2,m_l,generPoints,count=0,multiGauss=1):
	part =[0,0,0,0,0]
	sp = (m_Lb+m_Lc)**2-q2
	sm = (m_Lb-m_Lc)**2-q2
	formFack = all_f(q2,generPoints,count,multiGauss)
	part[0] = ((Gf**2 * sqrt(sp*sm))/(768*pi**3*m_Lb**3))*(1 - m_l**2/q2)**2
	part[1] = 4 * (m_l**2+2*q2) * (sp*((1-ei)*formFack['gt'])**2+sm*((1+ei)*formFack['ft'])**2)
	part[2] = 2*((m_l**2+2*q2)/q2)*(sp*((m_Lb-m_Lc)*(1-ei)*formFack['gp'])**2+sm*((m_Lb+m_Lc)*(1+ei)*formFack['fp'])**2)
	part[3] = 6*m_l**2/q2*(sp*((m_Lb-m_Lc)*(1+ei)*formFack['f0'])**2+sm*((m_Lb+m_Lc)*(1-ei)*formFack['g0'])**2)
	dG = part[0] * (part[1]+part[2]+part[3])
	return dG

def makeG(generPoints,count,multiGauss=1):
	q2=[0.01*i for i in range(1,1103)]
	Gm,Gt = [],[]
	for iQ in q2:
		Gm.append(partWidth(iQ,m_m,generPoints,count,multiGauss))
		Gt.append(partWidth(iQ,m_t,generPoints,count,multiGauss))
	return q2,Gm, Gt

def makeHistos(q2,Gm,Gt,count=""):
    hTau = r.TH1D("h_theor_tau"+count,"h_True_tau"+count,4,-2,14)
    hMu = r.TH1D("h_theor_mu"+count,"h_True_mu"+count,4,-2,14)
    for i in xrange(len(Gm)):
        if q2[i]>3.1: hTau.Fill(q2[i],Gt[i])
        hMu.Fill(q2[i],Gm[i])
  #  hTau.Scale(1.0/hTau.Integral("width"))
  #  hMu.Scale(1.0/hMu.Integral("width"))
    return hTau, hMu

def plotManyTemplates(plot=1):
    paramSample = obtaineTheorieParameters()
    hTauList, hMuList = [],[]
    for iSampl in xrange(len(paramSample[paramSample.keys()[0]])):
        q2,Gm, Gt = makeG(paramSample,iSampl)
        hTau, hMu = makeHistos(q2,Gm,Gt,str(iSampl))
        hTauList.append(hTau)
        hMuList.append(hMu)
    if plot:
        c=r.TCanvas()
        c.Divide(1,2)
        c.cd(1)
        for iH in xrange(20):#len(hTauList)):
            hTauList[iH].SetMarkerStyle(5)
            hMuList[iH].SetMarkerStyle(5)
            c.cd(1)
            hTauList[iH].SetMarkerColor(1+iH)
            if iH>0:hTauList[iH].Draw("SAME hist p")
            else:hTauList[iH].Draw("hist p")
            c.cd(2)
            if iH>0:hMuList[iH].Draw("SAME hist p")
            else:hMuList[iH].Draw("hist p")
            hMuList[iH].SetMarkerColor(1+iH)
        
    if plot: return hTauList, hMuList,c
    else: return hTauList, hMuList

def makeDistrPerQ2():
    hTauList, hMuList = plotManyTemplates(0)
    hTauQ2,hMuQ2 = {},{}
    for iQ in xrange(hTauList[0].FindFirstBinAbove(),hTauList[0].FindLastBinAbove()+1):
        hTauQ2[str(hTauList[0].GetBinCenter(iQ))] = r.TH1F("hTau_q2= "+str(hTauList[0].GetBinCenter(iQ)),"hTau_q2= "+str(hTauList[0].GetBinCenter(iQ)),3000,0,0.3)
    for iQ in xrange(hMuList[0].FindFirstBinAbove(),hMuList[0].FindLastBinAbove()+1):
        hMuQ2[str(hMuList[0].GetBinCenter(iQ))] = r.TH1F("hMu_q2= "+str(hMuList[0].GetBinCenter(iQ)),"hMu_q2= "+str(hMuList[0].GetBinCenter(iQ)),3000,0,0.3)
        
    for iH in hTauList:
        for iQ in xrange(iH.FindFirstBinAbove(),iH.FindLastBinAbove()+1):
            hTauQ2[str(iH.GetBinCenter(iQ))].Fill(iH.GetBinContent(iQ))
            
    for iH in hMuList:
        for iQ in xrange(iH.FindFirstBinAbove(),iH.FindLastBinAbove()+1):
            hMuQ2[str(iH.GetBinCenter(iQ))].Fill(iH.GetBinContent(iQ))
            
    return hTauQ2,hMuQ2

def makeDistrWithParamVariation():
   q2,Gm, Gt = makeG(0,0,0)
   hTauMean,hMuMean = makeHistos(q2,Gm,Gt)
   hTauQ2,hMuQ2 = makeDistrPerQ2()
   hBoundries = [r.TH1D("h_theor_"+iP+"_"+str(iT),"h_True_"+iP+"_"+str(iT),4,-2,14) for iT in ("max","min") for iP in("tau","mu")]
   for iH in hMuQ2:
       print float(hMuQ2[iH].GetName()[-4:])
       hBoundries[2].SetBinContent(hBoundries[2].FindBin(float(hMuQ2[iH].GetName()[-4:])),(hMuQ2[iH].GetMean()+hMuQ2[iH].GetStdDev()))
       hBoundries[3].SetBinContent(hBoundries[3].FindBin(float(hMuQ2[iH].GetName()[-4:])),(hMuQ2[iH].GetMean()-hMuQ2[iH].GetStdDev()))
       
   for iH in hTauQ2:   
       hBoundries[0].SetBinContent(hBoundries[0].FindBin(float(hTauQ2[iH].GetName()[-4:])),(hTauQ2[iH].GetMean()+hTauQ2[iH].GetStdDev()))
       hBoundries[1].SetBinContent(hBoundries[1].FindBin(float(hTauQ2[iH].GetName()[-4:])),(hTauQ2[iH].GetMean()-hTauQ2[iH].GetStdDev()))
   c=r.TCanvas()
   c.Divide(1,2)
   c.cd(1)
   
   hBoundries[0].SetFillStyle(3001)
   hBoundries[0].SetFillColorAlpha(r.kBlue, 0.85)
   hBoundries[1].SetFillColorAlpha(r.kWhite, 0.35)
   hBoundries[0].Draw()
   hTauMean.SetLineColor(r.kRed)
   hTauMean.SetFillStyle(3001)
   hTauMean.SetFillColorAlpha(r.kBlue, 0.85)
   hTauMean.Draw("hist SAME")
   hBoundries[1].Draw("SAME")
   c.cd(2)
   hBoundries[2].SetFillStyle(3001)
   hBoundries[2].SetFillColorAlpha(r.kBlue, 0.85)
  
   hBoundries[2].Draw()
   hBoundries[3].SetFillColorAlpha(r.kWhite, 0.35)
   
   hMuMean.SetLineColor(r.kRed)
   hMuMean.SetFillStyle(3001)
   hMuMean.SetFillColorAlpha(r.kBlue, 0.85)
   
   hMuMean.Draw("hist SAME")
   
   hBoundries[3].Draw("SAME")
   
   return hMuMean,hTauMean,hBoundries,c
    
   