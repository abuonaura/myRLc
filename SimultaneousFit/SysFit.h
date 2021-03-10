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

//vector<string> ch2fit = {"Isolated","Kenriched"};
using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace HistFactory;


class SysFit
{

	public:
		SysFit();
		~SysFit() {;}

		//Set fit configuration
		void AllowBarlowBeaston(Bool_t v) {BBeast = v; }       //Turn on Barlow-Beaston fit
		void CorrectFF_GS(Bool_t v)       {FF_GS = v; } 	   //Turn on FF corrections Lc Ground State 
		void CorrectFF_HES(Bool_t v)      {FF_HES = v; } 	   //Turn on FF corrections Lc Higher Excited State 
		void CorrectFF_LES(Bool_t v)      {FF_LES = v; } 	   //Turn on FF corrections Lc Higher Excited State 
		void CorrectSweights(Bool_t v)	  {swcorr = v; }       //Trun on sweight corrections

		void SetMCcategory(TString MCcat) {MCcategory = MCcat;} //Set MC category to fit
		void SetFitType(TString ftype)    {FitType = ftype;} //Set MC type of fit to perform
		
		void DoSimultaneousFit();      //Perform a simultaneous fit on the 3 channels
        void FitIsolated();            //Only perform fit on isolated ch.
        void FitLcpipi();              //Only perform fit on Lcpipi ch.
        void FitKenriched();		   //Only perform fit on Kenriched ch.
		void SelectChannel2fit(vector<string> ch) {channel_names = ch;} //Perform a simultaneous fit of selected channels

		void SetStartParameter(string, string, Double_t, Double_t, Double_t); //Set mean, min, max start value for a sample in a channel
		void SetTemplateFileName(string ch, string name);      //Sets name of files from which templates are read
		void SetSampleNames(TString);						   //Sets names of template samples for the MCcategory

		//Get fit configuration
		Bool_t IsBarlowBeastonOn() {return BBeast;}   //Retrieve BBeast value
		Bool_t IsGSFFcorrected()   {return FF_GS;}    //Retrieve if Lc GS is FF corrected
		Bool_t IsHESFFcorrected()  {return FF_HES;}   //Retrieve if Lc HES is FF corrected
		Bool_t IsLESFFcorrected()  {return FF_LES;}   //Retrieve if Lc LES is FF corrected
		Bool_t IsSweightCorrOn()   {return swcorr;}   //Retrieve if the sweight correction is applied 

		TString GetMCcategory() {return MCcategory;} //Retrieve which MC category is being fitted
		TString GetFitType(){return FitType;}		   //Retrieve Type of fit to perfom (single,isolated)

		map<string,vector<Double_t>> GetStartParameters (string); //Gets the starting Parameters for all samples in the chosen channel
		vector<Double_t> GetStartParameter (string, string);      //Gets the starting Parameters for one sample in the chosen channel

		string GetTemplateFileName(string ch);       		     //Gets name of files from which templates are read
		vector<string> GetSampleNames() {return sample_names;}   //Gets names of template samples considered
		vector<string> GetChannels2fit() {return channel_names;}    //return Isolated/Kenriched/Lcpipi
	

		//Functions for sweight correction
		TString GetSweightHistFileName()          {return swfname;}    			//Get the name of the file with templates sweight corrected
		void SetSweightHistFileName(TString fname){swfname = fname;}   			//Set the name of the file with templates sweight corrected
        void CreateSweightCorrectHistos(string, string, vector<string>); //Create the sweight corrected templates 
		//Functions for performing the fit
		
		Double_t GetHistoNormalisation(string, string);        		//Takes the normalization factor of each histogram to 1 (1./h->Integral())
		void AddSample(string,string, bool, const bool, Channel**); //Adds the samples to the channel
		ModelConfig* SetChannelConstants(ModelConfig*, string);          //Sets which parameters must be 
		HistFactory::Measurement CreateMeasurement();
																		 //constant for each channel

        RooWorkspace *CreateWorkspace(RooStats::HistFactory::Measurement);
        ModelConfig* CreateModel(RooWorkspace* );
        ModelConfig* FixYields(RooStats::ModelConfig*, string, string);
        ModelConfig* SetAlphaStartingPoints(RooStats::ModelConfig*);

		RooFitResult* Fit(ModelConfig*, Measurement, RooWorkspace *, Bool_t, Bool_t); //Perform the fit
		
		//Functions for getting fit results
		
		ModelConfig* LoadFitResults(string fname, ModelConfig* mc, string channelName);
        void SaveFitResults(string,RooFitResult *fitResult);
        //void StoreFitResults(string,RooFitResult *fitResult);
        //void CheckDiscrepancyWrtLastRLcValue(string fname, string chName);

		//Functions for plotting fit results
		
		//Gets name of the components for writing name in the legend
		TString GetComponentName(TString); 
		//Assigns each component entering the fit a colour
        Int_t GetComponentColor(TString); 
		//Returns the name of the fit variable to store the plot
		TString GetFitVarName(TString);
		//Plot all the three fit variables in one canvas
        void PlotFitVariables(RooRealVar* fitvar1,const char* title1, RooRealVar* fitvar2, const char* title2, RooRealVar* fitvar3, const char* title3,std::map<std::string,RooDataHist*> data ,RooSimultaneous* model, RooCategory* idx,string name_suffix);
		/*
		//Plot one fit variable
		void PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE);
		*/
        RooPlot* AdjustVarPlot(RooPlot *frame, double ymin_1, double ymax_1, double ymin_2, double ymax_2);
        RooPlot* AdjustPullsPlot(RooPlot* pframe, RooPlot *frame, double ymin_1, double ymax_1, double ymin_2, double ymax_2);
		/*
		//Plot one fit variable in bins of another variable
        void PlotInBins(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::ModelConfig *mc,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE);
*/


		//

	private:
		Bool_t BBeast;  //Turn on/off Barlow Beaston fit
		Bool_t FF_GS;   //Turn on/off FF corrections Lc Ground State (Lcmunu, Lctaunu)
		Bool_t FF_HES;  //Turn on/off FF corrections Lc Higher Excited State (Lc2625munu, Lc2625taunu)
		Bool_t FF_LES;  //Turn on/off FF corrections Lc Lower Excited State (Lc2593munu, Lc2593taunu)
		Bool_t swcorr;  //Turn on/off sweight corrections
		
		TString MCcategory;  //either MCfull or MCTrackerOnly
        TString FitType;     //Single, Simultaneous
		TString swfname;	 //Name file to read when applying sweights

		vector<string> channel_names; 			 //Isolated, Kenriched, Lcpipi
		map<string, string> template_fname;	     //Name file to read templates from
		vector<string> sample_names; 			 //Lb_Lcmunu, ...

		map<string,map<string,vector<Double_t>>> start_parameters; //Starting parameters for fit per category,sample,{value,min,max}

		RooStats::HistFactory::Measurement measure;
        RooWorkspace *wspace;
        RooStats::ModelConfig* model;
		
		//Parameters for blinding + blinding function
		Int_t alpha;
        Int_t beta ;
        Int_t gamma;
        Int_t alpha_s;
        Int_t beta_s;
        Int_t gamma_s;
		void blindResult(RooFitResult*,string name_suffix);



		

};


#endif
