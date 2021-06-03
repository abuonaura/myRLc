import ROOT as r
from ROOT import TMath as mt
import numpy as np

#Starting variables
m_pi = 139.57018 #+/- 0.00035 MeV (PDG)
m_K = 493.677 #+/- 0.016 MeV (PDG)
m_p = 938.272081 #+/- 0.000006 MeV (PDG)
m_mu = 105.6583745 #+/- 0.0000024 MeV (PDG)
m_Lc = 2286.46 #+/- 0.14 MeV (PDG)
ISOBDTcut=0.35
ISOBDT2cut=0.2
Kcodes = [130,310,311,321,313,323,10311,10321,10313,10323,20313,20323,30313,30323,10315,10325,20315,20325,317,327,319,329]
Dcodes = [411,421,10411,10421,413,423,10413,10423,20413,20423,415,425,431,10431,433,10433,20433,435]
picodes = [111,211]
pcodes = [2212]

def ScaleHisto(h,value):
    scale = value/h.Integral()
    h.Scale(scale)
    return h

def FillMotherHisto(h,pdg,w):
    if mt.Abs(pdg)==113: #rho0
        h.Fill(1,w)
    elif mt.Abs(pdg)==213: #rho+
        h.Fill(2,w)
    elif mt.Abs(pdg)==221: #eta
        h.Fill(3,w)
    elif mt.Abs(pdg)==331: #eta'
        h.Fill(4,w)
    elif mt.Abs(pdg)==223: #omega(782)
        h.Fill(5,w)
    elif mt.Abs(pdg)==333: #phi
        h.Fill(6,w)
    elif mt.Abs(pdg)==310: #K0S
        h.Fill(7,w)
    elif mt.Abs(pdg)==411: #D+
        h.Fill(8,w)
    elif mt.Abs(pdg)==421: #D0
        h.Fill(9,w)
    elif mt.Abs(pdg)==413: #D*+
        h.Fill(10,w)
    elif mt.Abs(pdg)==423: #D*0
        h.Fill(11,w)
    elif mt.Abs(pdg)==431: #Ds
        h.Fill(12,w)
    elif mt.Abs(pdg)==433: #D*s
        h.Fill(13,w)
    elif mt.Abs(pdg)==10433: #Ds1(2536)+
        h.Fill(14,w)
    elif mt.Abs(pdg)==20433:#Ds1(2460)+
        h.Fill(15,w)
    elif mt.Abs(pdg)==4122:#Lc
        h.Fill(16,w)
    elif mt.Abs(pdg)==104124:#Lc(2625)
        h.Fill(17,w)
    elif mt.Abs(pdg)==14122: #Lc(2593)
        h.Fill(18,w)
    elif mt.Abs(pdg)==5122: #Lb
        h.Fill(19,w)
    else:
        h.Fill(20,w)
    h.SetDirectory(0)
    return h

def CreateMotherHisto(name,title):
    h = r.TH1F(name,title, 20,1,21)
    MotherParticles = ['#rho^{0}','#rho^+','#eta','#eta^{1}','#omega(782)','#phi','K^{0}_{S}','D^{+}','D^{0}','D^{*+}','D^{*0}','D_{s}','D^{*}_{s}','D_{s1}(2536)^{+}','D_{s1}(2460)^{+}','#Lambda_{c}','#Lambda_{c}(2625)','#Lambda_{c}(2593)','#Lambda_{b}','others']
    for n,mp in enumerate(MotherParticles):
        h.GetXaxis().SetBinLabel(n+1,mp)
        h.GetXaxis().LabelsOption('v')
    h.SetDirectory(0)
    return h

def ComputeLcpipiMass(Lb_ISOLATION_PX, Lb_ISOLATION_PY, Lb_ISOLATION_PZ, Lb_ISOLATION_PX2, Lb_ISOLATION_PY2,Lb_ISOLATION_PZ2, Lc_PX, Lc_PY, Lx_PZ):
    #-> I assume that the 2 anti-isolated particles are 2 pions
    E1pi = mt.Sqrt(mt.Power(Lb_ISOLATION_PX,2)+mt.Power(Lb_ISOLATION_PY,2)+mt.Power(Lb_ISOLATION_PZ,2)+m_pi*m_pi)
    E2pi = mt.Sqrt(mt.Power(Lb_ISOLATION_PX2,2)+mt.Power(Lb_ISOLATION_PY2,2)+mt.Power(Lb_ISOLATION_PZ2,2)+m_pi*m_pi)
    ELc = mt.Sqrt(mt.Power(Lc_PX,2)+mt.Power(Lc_PY,2)+mt.Power(Lc_PZ,2)+m_Lc*m_Lc)
    #-> I compute the tot. energy of the 3 particle system
    pLc12_x = Lc_PX+Lb_ISOLATION_PX+Lb_ISOLATION_PX2
    pLc12_y = Lc_PY+Lb_ISOLATION_PY+Lb_ISOLATION_PY2
    pLc12_z = Lc_PZ+Lb_ISOLATION_PZ+Lb_ISOLATION_PZ2
    mLc12 = mt.Sqrt(mt.Power(ELc+E1pi+E2pi,2)-(mt.Power(pLc12_x,2)+mt.Power(pLc12_y,2)+mt.Power(pLc12_z,2)))
    return mLc12

def AssignMassIsoParticle_wTrueID(TrueID):
    m = m_pi
    isK=0
    if mt.Abs(TrueID) in Kcodes:
        m = m_K
        isK=1
    else:
        if mt.Abs(TrueID)==2212:
            m = m_p
        else:
            m = m_pi    
    return m, isK

def AssignMassIsoParticle_wPID(PIDK,muCh,Ch, PIPp):
    m = m_pi
    isK=0
    if PIDK>4.:
        if ((Ch==-muCh) or (Ch==muCh and (PIDp-PIDK)<=0.)):
            m = m_K
            isK=1
        if Ch==muCh and (PIDp-PIDK)>0.:
            m = m_p
    return m, isK

def ComputeLbMass(Lb_ISOLATION_PX, Lb_ISOLATION_PY, Lb_ISOLATION_PZ, Lb_ISOLATION_PX2, Lb_ISOLATION_PY2,Lb_ISOLATION_PZ2, Lc_PX, Lc_PY, Lx_PZ, TrueID, TrueID2, mu_PX, mu_PY, muPZ):
    m1, isK1 = AssignMassIsoParticle(TrueID)
    m2, isK2 - AssignMassIsoParticle(TrueID2)
    pLc12_x = Lc_PX+Lb_ISOLATION_PX+Lb_ISOLATION_PX2
    pLc12_y = Lc_PY+Lb_ISOLATION_PY+Lb_ISOLATION_PY2
    pLc12_z = Lc_PZ+Lb_ISOLATION_PZ+Lb_ISOLATION_PZ2
    E1 = mt.Sqrt(mt.Power(Lb_ISOLATION_PX,2)+mt.Power(Lb_ISOLATION_PY,2)+mt.Power(Lb_ISOLATION_PZ,2)+mt.Power(m1,2))
    E2 = mt.Sqrt(mt.Power(Lb_ISOLATION_PX2,2)+mt.Power(Lb_ISOLATION_PY2,2)+mt.Power(Lb_ISOLATION_PZ2,2)+mt.Power(m2,2))
    Emu = mt.Sqrt(mt.Power(mu_PX,2)+mt.Power(mu_PY,2)+mt.Power(mu_PZ,2)+m_mu*m_mu)
    ELc = mt.Sqrt(mt.Power(Lc_PX,2)+mt.Power(Lc_PY,2)+mt.Power(Lc_PZ,2)+m_Lc*m_Lc)
    ETOT = E1 + E2 + Emu + ELc
    pTOT_x = pLc12_x+t.mu_PX
    pTOT_y = pLc12_y+t.mu_PY
    pTOT_z = pLc12_z+t.mu_PZ
    mTOT = mt.Sqrt(mt.Power(ETOT,2)-pTOT_x*pTOT_x - pTOT_y*pTOT_y - pTOT_z*pTOT_z)
    return isK1, isK2, mTOT

def isKenriched(isK1,isK2,mTOT,mLc12):
    isKenr=0
    if mLc12>2700 and mTOT<5620 and (isK1==1 or isK2==1):
        isKenr=1
    return isKenr

def DiscardGhostsAndShortTracks(NNghost,NNghost2,Type,Type2):
    accept = 0
    if Type==3 and NNghost<0.2 and Type2==3 and NNghost2<0.2:
        accept=1
    return accept

def PlotHisto(h,name):
    c = r.TCanvas(name,'',500,500)
    h.SetLineColor(r.kAzure-3)
    h.Draw('hist')
    return c

def PlotHistoSameCanvas(c,h,color):
    h.SetLineColor(color)
    h.Draw('hist same')
    return c

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


def SelectTrueKevents(TrueID, TrueMID, TrueID2, TrueMID2):
    Kfound, Kfound2=0,0
    KfromDfound,KfromDfound2=0,0
    if mt.Abs(TrueID) in Kcodes:
        Kfound=1
        if mt.Abs(TrueMID) in Dcodes:
            KfromDfound=1
    if mt.Abs(TrueID2) in Kcodes:
        Kfound2=1
        if mt.Abs(TrueMID2) in Dcodes:
            KfromDfound2=1
    if Kfound==1 or Kfound2==1:
        if KfromDfound==1 or KfromDfound2==1:
            return 1
        else:
            return 0
    else:
        return 0
            
def SelectKfromDevents(TrueID, TrueMID):
    if mt.Abs(TrueID) in Kcodes:
        if mt.Abs(TrueMID) in Dcodes:
            return 1
        else:
            return 0
    else:
        return 0

#Upload files to read
datadir = '/disk/lhcb_data2/RLcMuonic2016/MC_full_new/'
#samples = ['Lcmunu','Lctaunu','LcDs','Lc2593munu','Lc2593taunu','Lc2593Ds','Lc2625munu','Lc2625taunu','Lc2625Ds']
samples = ['LcDs']

for sample in samples:
    print('--------------- '+sample+' ---------------')
    f = r.TFile(datadir+'Lb_'+sample+'_MagUp.root','READ')
    f1 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_preselectionVars.root')
    f2 = r.TFile(datadir+'Lb_'+sample+'_MagUp_full_PIDCalib.root')

    #Get Tree + tree with preselection variables
    t = f.Get('tupleout/DecayTree')
    t1 = f1.Get('DecayTree')
    t2 = f2.Get('tupleout/DecayTree')
    t.AddFriend(t1)
    t.AddFriend(t2)
    #Switch off all branches to speed up
    t.SetBranchStatus('*',0)
    #List of branches to turn on
    blist = ['FinalSel','isKenriched','isIsolated','Lb_ISOLATION_TruePID','Lb_ISOLATION_TrueMotherPID','Lb_ISOLATION_TruePID2','Lb_ISOLATION_TrueMotherPID2','Lb_ISOLATION_BDT','Lb_ISOLATION_BDT2','Lb_ISOLATION_BDT3','Lb_ISOLATION_PIDK','Lb_ISOLATION_PIDp','Lb_ISOLATION_PIDK2','Lb_ISOLATION_PIDp2','Lb_ISOLATION_PIDK3','Lb_ISOLATION_PIDp3','Lb_ISOLATION_Type', 'Lb_ISOLATION_Type2','Lb_ISOLATION_Type3','Lb_ISOLATION_NNghost', 'Lb_ISOLATION_NNghost2','Lb_ISOLATION_NNghost3', 'Lb_ISOLATION_CHARGE', 'Lb_ISOLATION_CHARGE2','Lb_ISOLATION_CHARGE3','Event_PIDCalibEffWeight','mu_ID','Lb_ISOLATION_PX','Lb_ISOLATION_PX2','Lb_ISOLATION_PY','Lb_ISOLATION_PY2','Lb_ISOLATION_PZ','Lb_ISOLATION_PZ2','Lc_PX','Lc_PY','Lc_PZ','mu_PX','mu_PY','mu_PZ','FitVar_q2_mLc','FitVar_Mmiss2_mLc','FitVar_El_mLc']
    for b in blist:
        t.SetBranchStatus(b,1)

    hBDTall = r.TH1F('hBDTall','ISOBDT',25,-1,1) #ISOBDT all anti-iso particle
    hBDT2all = r.TH1F('hBDT2all','ISOBDT2',25,-1,1)
    hMotherIDall = CreateMotherHisto('hMotherIDall','Mother ID all anti-iso particles')
    hMotherID2all = CreateMotherHisto('hMotherID2all','Mother ID all anti-iso particles (2)')
    hMotherIDK = CreateMotherHisto('hMotherIDK','Mother ID True K')
    hMotherID2K = CreateMotherHisto('hMotherID2K','Mother ID True K (2)')
    hBDTKfromD = r.TH1F('hBDTKfromD','ISOBDT of K from D',25,-1,1) #ISOBDT true K from charm hadrons
    hBDTKfromD2 = r.TH1F('hBDTKfromD2','ISOBDT2 of K from D',25,-1,1) #ISOBDT2 true K from charm hadrons
    h_templ_allK    = r.TH3F('h_templ_allK',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)
    h_templ_KfromD = r.TH3F('h_templ_KfromD',';q^{2};E_{l};M^{2}_{miss}',4,-2,14,10,0,2600,10,-2,14)


    #Loop over entries
    #for i in range(t.GetEntries()):
    for i in range(500000):
        t.GetEntry(i)
        #discard events not passing all the selection requirements
        if t.FinalSel!=1:
            continue
        muCh = -t.mu_ID/13
        BDT = t.Lb_ISOLATION_BDT
        PIDp = t.Lb_ISOLATION_PIDp
        PIDK = t.Lb_ISOLATION_PIDK
        Ch   = t.Lb_ISOLATION_CHARGE
        Type = t.Lb_ISOLATION_Type
        TrueID = t.Lb_ISOLATION_TruePID
        TrueMID = t.Lb_ISOLATION_TrueMotherPID
        NNghost = t.Lb_ISOLATION_NNghost
        BDT2 = t.Lb_ISOLATION_BDT2
        PIDp2 = t.Lb_ISOLATION_PIDp2
        PIDK2 = t.Lb_ISOLATION_PIDK2
        Ch2 = t.Lb_ISOLATION_CHARGE2
        Type2 = t.Lb_ISOLATION_Type2
        TrueID2 = t.Lb_ISOLATION_TruePID2
        TrueMID2 = t.Lb_ISOLATION_TrueMotherPID2
        NNghost2 = t.Lb_ISOLATION_NNghost2
        w = t.Event_PIDCalibEffWeight #PID calib reweight

        #Discard ghost tracks (NNghost>=0.2 and not long tracks !=3)
        if(DiscardGhostsAndShortTracks(NNghost,NNghost2,Type,Type2)==1):
            hBDTall.Fill(BDT,w)
            hBDT2all.Fill(BDT2,w)
            #print(TrueMID)
            hMotherIDall = FillMotherHisto(hMotherIDall,TrueMID,w)
            hMotherID2all = FillMotherHisto(hMotherID2all,TrueMID2,w)
            if TrueID in Kcodes:
                hMotherIDK = FillMotherHisto(hMotherIDK,TrueMID,w)
                h_templ_allK.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
            if TrueID2 in Kcodes:
                hMotherID2K = FillMotherHisto(hMotherID2K,TrueMID2,w)
                h_templ_allK.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
            #Select K coming from D decays:
            #if SelectTrueKevents(TrueID, TrueMID, TrueID2, TrueMID2):
            if SelectKfromDevents(TrueID, TrueMID):
                hBDTKfromD.Fill(BDT,w)
                h_templ_KfromD.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)
            if SelectKfromDevents(TrueID2, TrueMID2):
                hBDTKfromD2.Fill(BDT2,w)
                h_templ_KfromD.Fill(t.FitVar_q2_mLc/1.E6,t.FitVar_El_mLc,t.FitVar_Mmiss2_mLc/1.E6,w)



            
    c = PlotHisto(hMotherIDall,'c')
    c = PlotHistoSameCanvas(c, hMotherIDK,r.kOrange-3)
    c1 = PlotHisto(hMotherID2all,'c1')
    c1 = PlotHistoSameCanvas(c1, hMotherID2K,r.kOrange-3)
    c2 = PlotHisto(hBDTall,'c2')
    c2 = PlotHistoSameCanvas(c2,hBDTKfromD,r.kOrange-3)
    c3 = PlotHisto(hBDT2all,'c3')
    c3 = PlotHistoSameCanvas(c3,hBDTKfromD2,r.kOrange-3)
    c4 = PlotTemplates(h_templ_allK,'c4')
    c4 = PlotTemplatesSameCanvas(c4,h_templ_KfromD,r.kOrange-3)
