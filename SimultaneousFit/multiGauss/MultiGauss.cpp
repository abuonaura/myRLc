#include "MultiGauss.h"
#include <iostream>
namespace multifit {

    MultiGauss::MultiGauss(){}
    
    MultiGauss::~MultiGauss(){}
    
    map<Int_t, std::string> MultiGauss::stupidNaming()
    {
        map<Int_t, std::string> stupidMap;
        stupidMap[0] = {"a0fp"};
        stupidMap[1] = {"a1fp"};
        stupidMap[2] = {"a0f0"};
        stupidMap[3] = {"a1f0"};
        stupidMap[4] = {"a0ft"};
        stupidMap[5] = {"a1ft"};
        stupidMap[6] = {"a0gp"};
        stupidMap[7] = {"a1gp"};
        stupidMap[8] = {"a0g0"};
        stupidMap[9] = {"a1g0"};
        stupidMap[10] = {"a1gt"};
        return stupidMap;
    }
    unordered_map <std::string, std::vector<Double_t> > MultiGauss::defineValues()
    {
        unordered_map<std::string, std::vector<Double_t> > parVal;
        
        parVal["a0fp"] = {0.8146,0.0167};
        parVal["a1fp"] = {-4.8990,0.5425};
        parVal["a0f0"] = {0.7439,0.0125};
        parVal["a1f0"] = {-4.6480,0.6084};
        parVal["a0ft"] = {1.0780,0.0256};
        parVal["a1ft"] = {-6.4170,0.8480};
        parVal["a0gp"] = {0.6847,0.0086};
        //parVal["a0gt"] = {0.6847,0.0086};
        parVal["a1gp"] = {-4.4310,0.3572};
        parVal["a0g0"] = {0.7396,0.0143};
        parVal["a1g0"] = {-4.3660,0.3314};
        parVal["a1gt"] = {-4.4630,0.3613};
        return parVal;
    }

    TMatrixDSym MultiGauss::covMat(map<Int_t, std::string> stupidMap,unordered_map <std::string, std::vector<Double_t> > values)
    {   
        Int_t dim = 11;
        TMatrixDSym corr(dim);
        for(Int_t i=0;i<11;i++){corr(i,i)=1;}
            corr(0,1) = corr(1,0) = -0.6644;
            corr(0,2) = corr(2,0) = 0.6827;
            corr(1,2) = corr(2,1) = -0.6515;
            corr(0,3) = corr(3,0) = -0.4853;
            corr(1,3) = corr(3,1) = 0.9445;
            corr(2,3) = corr(3,2) = -0.7040;
            
            corr(0,4) = corr(4,0) = 0.6218;
            corr(1,4) = corr(4,1) = -0.3853;
            corr(2,4) = corr(4,2) = 0.4208;
            corr(3,4) = corr(4,3) = -0.2738;
            
            corr(0,5) = corr(5,0) = -0.3906;
            corr(1,5) = corr(5,1) = 0.5109;
            corr(2,5) = corr(5,2) = -0.3620;
            corr(3,5) = corr(5,3) = 0.4739;
            corr(4,5) = corr(5,4) = -0.6637;
 
            corr(0,6) = corr(6,0) = 0.4828;
            corr(1,6) = corr(6,1) = -0.3831;
            corr(2,6) = corr(6,2) = 0.6174;
            corr(3,6) = corr(6,3) = -0.3888;
            corr(4,6) = corr(6,4) = 0.3933;
            corr(5,6) = corr(6,5) = -0.2903;
 
            corr(0,7) = corr(7,0) = -0.3152;
            corr(1,7) = corr(7,1) = 0.4915;
            corr(2,7) = corr(7,2) = -0.4822;
            corr(3,7) = corr(7,3) = 0.5261;
            corr(4,7) = corr(7,4) = -0.2369;
            corr(5,7) = corr(7,5) = 0.3509;
            corr(6,7) = corr(7,6) = -0.7304;
 
            corr(0,8) = corr(8,0) = 0.5636;
            corr(1,8) = corr(8,1) = -0.2979;
            corr(2,8) = corr(8,2) = 0.4320;
            corr(3,8) = corr(8,3) = -0.2164;
            corr(4,8) = corr(8,4) = 0.5161;
            corr(5,8) = corr(8,5) = -0.2443;
            corr(6,8) = corr(8,6) = 0.6365;
            corr(7,8) = corr(8,7) = -0.3829;
         
            
            corr(0,9) = corr(9,0) = -0.4317;
            corr(1,9) = corr(9,1) = 0.4916;
            corr(2,9) = corr(9,2) = -0.4726;
            corr(3,9) = corr(9,3) = 0.4779;
            corr(4,9) = corr(9,4) = -0.3639;
            corr(5,9) = corr(9,5) = 0.3640;
            corr(6,9) = corr(9,6) = -0.6743;
            corr(7,9) = corr(9,7) = 0.8725;
            corr(8,9) = corr(9,8) = -0.6843;
            
            
            corr(0,10) = corr(10,0) = -0.3763;
            corr(1,10) = corr(10,1) = 0.4764;
            corr(2,10) = corr(10,2) = -0.4756;
            corr(3,10) = corr(10,3) = 0.4877;
            corr(4,10) = corr(10,4) = -0.2926;
            corr(5,10) = corr(10,5) = 0.3400;
            corr(6,10) = corr(10,6) = -0.7301;
            corr(7,10) = corr(10,7) = 0.9171;
            corr(8,10) = corr(10,8) = -0.4846;
            corr(9,10) = corr(10,9) = 0.8456;
            
            TMatrixDSym cov(dim);
            for(Int_t i=0;i<dim;i++){
                for(Int_t j=0;j<dim;j++){
                    cov(i,j)=corr(i,j)*values[stupidMap[i]][1]*values[stupidMap[j]][1];
                }
            }
        return cov;
    }
    map<std::string, std::vector<Double_t> >  MultiGauss::function()
    {
        RooArgList params,unc;
        RooArgSet poi;
        unordered_map<std::string, std::vector<Double_t> > parVal = defineValues();
        map<Int_t, std::string> stupidMap = stupidNaming();
        TMatrixDSym cov = covMat(stupidMap,parVal);
        for (map<Int_t, std::string>::iterator iP= stupidMap.begin(); iP != stupidMap.end();++iP){
              params.add(*new RooRealVar(Form("x%s",iP->second.c_str()),Form("x%s",iP->second.c_str()),0,-20,20));
              unc.add(*new RooRealVar(Form("x_un%s",iP->second.c_str()),Form("x_un%s",iP->second.c_str()),parVal[iP->second][0]));
        }
        
        RooMultiVarGaussian mvg("mvg", "mvg", params, unc, cov);
        RooDataSet* data = mvg.generate(params, 1000);
        map<std::string, std::vector<Double_t> > dataPoints;
        TH1F* h =new TH1F("h","h",200,-10,10);
        for(Int_t i=0;i<1000;i++){
            const RooArgSet* row = data->get(i);
            for (unordered_map<std::string, std::vector<Double_t> >::iterator iP= parVal.begin(); iP != parVal.end();++iP){
            RooRealVar* xrow = (RooRealVar*) row->find(Form("x%s",iP->first.c_str()));
            dataPoints[Form("x%s",iP->first.c_str())].push_back(xrow->getVal());   
            }
            RooRealVar* temp = (RooRealVar*) row->find("xa0f0");
            h->Fill(temp->getVal());
        }
        h->Draw();
        return  dataPoints;
    }

}
