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


class SysFit{

	public:
		SysFit();
		~SysFit() {;}
		vector<string> NameChannels(); //return Isolated/Kenriched
		void AllowBarlowBeaston() {BBeast = true;} //if this function is called the BarlowBeaston fit is turned on
		map<string,vector<Double_t>> GetStartParameters (string); //Gets the starting Parameters for the different samples in the 2 channels 
		void SetStartParameters(map<string,vector<Double_t> >); //sets the start fit parameters
		void PrintStartParams(string, map<string,vector<Double_t>>); //prints the starting fit parameters
		vector<string> GetCategory(map<string,vector<Double_t> >); //Retrieves the name of the samples + "pha" (?)
		TString GetComponentName(TString); //Function to call when creating the legend
        	Int_t GetComponentColor(TString); //Assigns to each sample a color
		TString GetFitVarName(TString); 
		Double_t GetHistoNormalisation(string, string); //Takes the normalization factor of each histogram to 1 (1./h->Integral())
		Bool_t ActivateShapeUncertainties(string,RooStats::HistFactory::Channel*);
		Bool_t ActivateGaussianConstraint(string,RooStats::HistFactory::Channel*);
		void AddSample(string,string, bool ,const bool, bool, RooStats::HistFactory::Channel**,vector<Double_t>,vector<Double_t>); //Adds the samples to the channel
		vector<Double_t> Fit(); //Perform the fit
		void PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::HistFactory::HistFactorySimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, bool legend=kFALSE);


	private:
		//yields for applying gaussian constraints for MISID/Combinatorial in the control fit
        Double_t y_misid_cf;
        Double_t y_comb_cf ;
        //weights for applying gaussian constraints for MISID/Combinatorial in the control fit
        Double_t w_misid_cf;
        Double_t w_comb_cf;

        Int_t alpha;
        Int_t beta ;
        Int_t gamma;
        Int_t alpha_s;
        Int_t beta_s;
        Int_t gamma_s;

	Bool_t BBeast;
	Bool_t ShapeUnc;
	Bool_t GaussConstr;

	map<string,vector<Double_t>> start_parameters;
	map<string,vector<Double_t> > fit_result;


};


#endif
