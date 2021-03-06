'''
Author: Iaroslava Bezshyiko
Modified by: Annarita Buonaura

'''

import ROOT as r
import math as m
import array
import argparse


def getDistr():
    coefFile = r.TFile.Open("PTETAweights_2016_Aug18_250b_L0M.root")
    hCoefficients = coefFile.Get("hFold0")
    return hCoefficients

def countPseudorapidity(t):
    P_TRUE = m.sqrt(t.Lb_TRUEP_X*t.Lb_TRUEP_X+t.Lb_TRUEP_Y*t.Lb_TRUEP_Y+t.Lb_TRUEP_Z*t.Lb_TRUEP_Z)
    if(P_TRUE!=0):
        eta = m.atanh(float(t.Lb_TRUEP_Z)/float(P_TRUE))
    else:
        eta = -99
    return eta
    
def correctLb(fileName,treeName):
    coefFile = r.TFile.Open("PTETAweights_2016_Aug18_250b_L0M.root")
    hCoefficients = coefFile.Get("hFold0")
    fileCorr = r.TFile.Open(fileName)
    treeCorr = fileCorr.Get(treeName)
    fileNew = r.TFile.Open(fileName[:-5]+"_LbCorrection.root","RECREATE")
    treeNew = r.TTree("DecayTree","DecayTree")
    w_LbCorr=array.array('f',1*[0])
    eta=array.array('f',1*[0])
    treeNew.Branch('w_LbCorr',w_LbCorr,'w_LbCorr[1]/F')
    treeNew.Branch('eta',eta,'eta[1]/F')
    i=0
    print(treeCorr.GetEntries())
    for iE in treeCorr:
        i+=1
        if i%50000==0: print( "  Event: ", i)
        eta[0] = countPseudorapidity(iE)
        if eta[0]==-99:
            w_LbCorr[0]=0
        else:
            pt = iE.Lb_TRUEPT*(0.0001)
            if pt>50000:
                pt = 49999*0.0001
            w_LbCorr[0] = hCoefficients.GetBinContent(hCoefficients.FindBin(eta[0],pt))
        #print(eta[0], w_LbCorr[0])
        treeNew.Fill()
    print(treeNew.GetEntries(), i)
    fileNew.Write()
    fileNew.Close()
    fileCorr.Close()
    coefFile.Close()
    


    
def compareReweight():
   # iFL = r.TFile.Open("DemoHistosLattice.root")
   # iF = r.TFile.Open("DemoHistos.root")
    var = {'P':['Lb_P',[40,0,6e5]],'Pt':['Lb_PT',[100,0,4e4]],'mu_Pt':['mu_PT',[40,0,3e4]],"IPCH2":["Lb_IPCHI2_OWNPV",[100,0,5e3]],"vtx2":["Lb_ENDVERTEX_CHI2",[100,0,6]]}
    types = {
    'mu':{'file':['Lb_Lcmunu_PID_reduced_preselected_LbCorr.root'],'tree':'DecayTree','weight':'Lc_BKGCAT==0&&Lb_BKGCAT<60','rate':'(0.584)'},
    'tau':{'file':['Lb_Lctaunu_PID_reduced_preselected_LbCorr.root'],'tree':'DecayTree','weight':'Lc_BKGCAT==0&&Lb_BKGCAT<60','rate':'(0.033)'},
    'Ds':{'file':['Lb_LcDs_PID_reduced_preselected_LbCorr.root'],'tree':'DecayTree','weight':'Lc_BKGCAT==0&&Lb_BKGCAT<60','rate':'(0.113)'},
    'LcStarmu':{'file':['Lb_LcStarmunu_PID_reduced_preselected_LbCorr.root','Lb_Lc2625munu_PID_reduced_preselected_LbCorr.root'],'tree':'DecayTree','weight':'Lc_BKGCAT==0&&Lb_BKGCAT<60','rate':'(0.188)'},
    'LcStartau':{'file':['Lb_LcStartaunu_PID_reduced_preselected_LbCorr.root'],'tree':'DecayTree','weight':'Lc_BKGCAT==0&&Lb_BKGCAT<60','rate':'(0.009)'},
    'data':{'file':['Lb_Data_MagUp_reduced_preselected_iso_sw.root','Lb_Data_MagDown_reduced_preselected_iso_sw.root'],'tree':'DecayTree','weight':'(sw_sig)'}
    }
    
    hMcDataNoCorr = {iV:{iT:r.TH1D("NodataMCnC"+iT+iV,"NodataMCnC"+iT+iV,var[iV][1][0],var[iV][1][1],var[iV][1][2]) for iT in types.keys()} for iV in var.keys() }
    hMcData = {iV:{iT:r.TH1D("dataMC"+iT+iV,"dataMC"+iT+iV,var[iV][1][0],var[iV][1][1],var[iV][1][2]) for iT in types.keys()} for iV in var.keys()}
    
    hData = {iV:r.TH1D("data"+iV,"data"+iV,var[iV][1][0],var[iV][1][1],var[iV][1][2]) for iV in var.keys()}
    hMC = {iV:r.TH1F('hMC_'+iV,'hMC_'+iV,var[iV][1][0],var[iV][1][1],var[iV][1][2]) for iV in var.keys()}
    hMCN = {iV:r.TH1F('hMCN'+iV,'hMCN'+iV,var[iV][1][0],var[iV][1][1],var[iV][1][2]) for iV in var.keys()}
    
    
    for iT in types.keys():
         if iT=='data': continue
         tree=r.TChain(types[iT]['tree'])
         for iF in types[iT]['file']:
            tree.Add(iF)
            
         for iV in var.keys():
          tree.Draw(var[iV][0]+">>dataMC"+iT+iV+"("+str(var[iV][1][0])+","+str(var[iV][1][1])+","+str(var[iV][1][2])+")","("+types[iT]['weight']+')'+'*'+types[iT]['rate']+"*w_LbCorr")
          hMcData[iV][iT] = r.gPad.GetPrimitive("dataMC"+iT+iV)
      #   hMcData[iT].Scale(1794766./hMcData[iT].Integral())
    for iV in hMcData.keys():
        for hi in hMcData[iV].keys():
            if hi=='data': continue
            hMC[iV].Add(hMcData[iV][hi])
        hMC[iV].Scale(1338263./hMC[iV].Integral())#1457842./hMC[iV].Integral())#1794766./hMC.Integral())
    
    for iT in types.keys():
         if iT=='data': continue
         tree=r.TChain(types[iT]['tree'])
         for iF in types[iT]['file']:
            tree.Add(iF)
         for iV in var.keys():
            tree.Draw(var[iV][0]+">>NodataMCnC"+iT+iV+"("+str(var[iV][1][0])+","+str(var[iV][1][1])+","+str(var[iV][1][2])+")","("+types[iT]['weight']+')*'+types[iT]['rate'])
            hMcDataNoCorr[iV][iT] = r.gPad.GetPrimitive("NodataMCnC"+iT+iV)
            hMcDataNoCorr[iV][iT].Scale(1./hMcDataNoCorr[iV][iT].Integral())
    for iV in hMcData.keys():
        for hi in hMcData[iV].keys():
            if hi=='data': continue
            hMCN[iV].Add(hMcDataNoCorr[iV][hi])
        hMCN[iV].Scale(1338263./hMCN[iV].Integral())#1457842./hMCN[iV].Integral())#1794766./hMCN.Integral())
    
    tree=r.TChain(types['data']['tree'])
    for i in types['data']['file']:
            tree.Add(i)
    for iV in var.keys():       
        tree.Draw(var[iV][0]+">>data_"+iV+"("+str(var[iV][1][0])+","+str(var[iV][1][1])+","+str(var[iV][1][2])+")",types['data']['weight'])
        hData[iV] = r.gPad.GetPrimitive("data_"+iV)
        hMC[iV].SetDirectory(0)
        hMCN[iV].SetDirectory(0)
        hData[iV].SetDirectory(0)
    return hMCN,hMC,hData
    
def draw():
    hMCN,hMC,hData = compareReweight()
    c=r.TCanvas()
    c.Divide(2,3)
    for counter,iV in enumerate(hData.keys()):
        print (counter)
        hMCN[iV].SetLineStyle(2)
        hData[iV].SetMarkerStyle(8)
        hMCN[iV].SetLineColor(r.kRed)
        hMC[iV].SetLineColor(r.kGreen)
        hMCN[iV].SetLineWidth(3)
        hMC[iV].SetLineWidth(3)
        hMC[iV].SetTitle(iV)
        c.cd(counter+1)
        hMC[iV].Draw("hist")
        hMCN[iV].Draw("same hist")
        hData[iV].Draw("same hist P")
    c.cd(6)
    legend=r.TLegend(0.1,0.7,0.48,0.9)
    legend.SetHeader("Lb correction","C")
    legend.AddEntry(hMCN["P"],"MC before correction","l");
    legend.AddEntry(hMC["P"],"MC after correction","l");
    legend.AddEntry(hData["P"],"Data","p");
    legend.Draw();
    return c, hMCN,hMC,hData,legend
    



