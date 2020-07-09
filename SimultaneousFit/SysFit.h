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

vector<string> ch2fit = {"Isolated","Kenriched"};
using namespace std;
using namespace RooFit;
using namespace RooStats;


class SysFit{

	public:
		SysFit();
		~SysFit() {;}

		void DoSimultaneousFit();
		void FitIsolated();
		void FitKenriched();
		string GetFitType(){return FitType;}
		void SetMCcathegory(string MCcat) {MCcathegory = MCcat;}
		string GetMCcathegory() {return MCcathegory;}
		void SelectChannel2fit(vector<string> ch) {channel_names = ch;}
		vector<string> NameChannels() {return channel_names;} //return Isolated/Kenriched
		
		void AllowBarlowBeaston() {BBeast = true;} //if this function is called the BarlowBeaston fit is turned on
		void ActivateShapeUncertainties(string, Bool_t);
		Bool_t IsShapeUncertain(string);
		void ActivateGaussConstraint(string, Bool_t);
		Bool_t IsGaussConstrained(string);
		void SetWeightGaussConstraint(string , Double_t );
		Double_t GetWeightGaussConstraint(string);

		map<string,vector<Double_t>> GetStartParameters (string); //Gets the starting Parameters for the different samples in the 2 channels 
		void SetStartParameters(RooFitResult *, string); //sets the start fit parameters
		void PrintStartParams(string, map<string,vector<Double_t>>); //prints the starting fit parameters

		vector<string> GetCategory(map<string,vector<Double_t> >); //Retrieves the name of the samples + "pha" (?)
		vector<string> GetParametersName(map<string,vector<Double_t> >); //Retrieves the name of the parameters + "pha" (?)
		
		TString GetComponentName(TString); //Function to call when creating the legend
        	Int_t GetComponentColor(TString); //Assigns to each sample a color
		TString GetFitVarName(TString); 
		
		Double_t GetHistoNormalisation(string, string); //Takes the normalization factor of each histogram to 1 (1./h->Integral())

		RooStats::ModelConfig* SetChannelConstants(RooStats::ModelConfig*, string);
		
		void AddSample(string,string, bool ,const bool, bool, RooStats::HistFactory::Channel**,vector<Double_t>,vector<Double_t>); //Adds the samples to the channel

		RooStats::HistFactory::Measurement CreateMeasurement();
		RooWorkspace *CreateWorkspace(RooStats::HistFactory::Measurement);
		RooStats::ModelConfig* CreateModel(RooWorkspace* );
		RooStats::ModelConfig* FixYields(RooStats::ModelConfig*, string, string);


		RooFitResult* Fit(RooStats::ModelConfig*, RooStats::HistFactory::Measurement, RooWorkspace *); //Perform the fit
		//void PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::HistFactory::HistFactorySimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE);
		void PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE);
		void PlotInBins(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::ModelConfig *mc,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE);

		void SaveFitResults(string,RooFitResult *fitResult);	
		void CheckDiscrepancyWrtLastRLcValue(string fname);

	private:

        Int_t alpha;
        Int_t beta ;
        Int_t gamma;
        Int_t alpha_s;
        Int_t beta_s;
        Int_t gamma_s;

	vector<string> channel_names;
	string MCcathegory; //either MCfull or MCTrackerOnly
	string FitType; //Single, Simultaneous

	Bool_t BBeast;
	

	map<string,map<string,vector<Double_t>>> start_parameters;
	map<string,vector<Double_t> > fit_result;
	map<string, Double_t> weight; //Weight for Gaussian Constraint 
	map<string, Bool_t> ShapeUnc; //set or not Gaussian Constraint 
	map<string, Bool_t> GaussConstr; //set or not Gaussian Constraint 

RooStats::HistFactory::Measurement measure;
	RooWorkspace *wspace;
	RooStats::ModelConfig* model;

	void blindResult(RooFitResult*);

};


#endif
