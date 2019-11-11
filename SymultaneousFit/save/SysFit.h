#ifndef SYSFIT_H
#define SYSFIT_H

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TLegend.h"

#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooMCStudy.h"
#include "RooMinuit.h"
#include "RooCategory.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooMsgService.h"
#include "RooParamHistFunc.h"
#include "RooHist.h"
#include "RooRandom.h"
#include "RooPlot.h"
#include "RooCmdArg.h"
#include "RooCustomizer.h"
#include "RooErrorVar.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/HistFactorySimultaneous.h"
#include "RooStats/HistFactory/Channel.h"
#include "RooStats/HistFactory/MakeModelAndMeasurementsFast.h"
#include "RooStats/HistFactory/Measurement.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooStats/HistFactory/HistFactoryModelUtils.h"
#include "RooStats/HistFactory/RooBarlowBeestonLL.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/MinNLLTestStat.h"
#include "RooStats/HistFactory/Sample.h"

#include <vector>
#include <map>


using namespace std;
using namespace RooFit;
using namespace RooStats;

//namespace RLC 
//{
    class SysFit
    {
    public:
        SysFit();
        ~SysFit();
        
	void SetNoOutputOnScreen(bool sw){
            if(sw){gROOT->SetBatch(kTRUE);}
            else{gROOT->SetBatch(kFALSE);}
        };
        
	vector<Double_t> Fit();
        Double_t GetHistoNormalisation(string, string);
        
	TString GetComponentName(TString);
	Int_t GetComponentColor(TString);
	TString GetFitVarName(TString);
	void PlotFrame(RooRealVar*,const char*,RooAbsData*,RooStats::HistFactory::HistFactorySimultaneous*,   RooCategory*,Double_t, Double_t, const char*, bool);
        void SetBB(bool bb){BBeast = bb;};
        std::map<std::string,std::vector<Double_t> > GetFitResult();
        map<string,vector<Double_t> > GetStartParameters(string);
        void SetStartParameters(std::map<std::string,std::vector<Double_t> >);
        void RLcPlot();
	vector<string> NameChannels();
	void PrintStartParams(string, map<string,vector<Double_t>>);
    
    private:
        bool BBeast = false;
        map<string,vector<Double_t> > fit_result;
        void AddSample(string,string, bool ,const bool, bool, RooStats::HistFactory::Channel*,vector<Double_t>,vector<Double_t>);
        RooFitResult *fitResult;
        Double_t correlation;
        void drawEllipse(Double_t,Double_t,Double_t,Double_t,Int_t, Double_t);
        void blindResult();
	
	//yields for applying gaussian constraints for MISID/Combinatorial in the control fit
	Double_t y_misid_cf = 5.E3;
	Double_t y_comb_cf = 7.E3;
	//weights for applying gaussian constraints for MISID/Combinatorial in the control fit
	Double_t w_misid_cf = 29.;
	Double_t w_comb_cf = 109.;
        
	Int_t alpha = 17;
        Int_t beta = 34;
        Int_t gamma = 98;
        Int_t alpha_s = 456;
        Int_t beta_s = 172;
        Int_t gamma_s = 23;
    };
//}

#endif
