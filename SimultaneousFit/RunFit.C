#include "SysFit.h"
#include "SysFit.cpp"
#include <TROOT.h>
#include <vector>
#include <iterator>
#include <map>
#include <iostream>

using namespace std;

Bool_t BBeast    =  0;
Bool_t FFGstate  =  1;
Bool_t FFEstateL =  0;
Bool_t FFEstateH =  1; 
Bool_t SWcorr    =  1; 
Bool_t MCTO      =  0; //Uses MC TrackerOnly samples
Bool_t FitIso    =  1; //Fit Isolated cat
Bool_t FitKenr   =  0; //Fit Kenriched cat
Bool_t FitLcpi   =  0; //Fit Lcpipi cat


void PrintFitConfiguration(string MCcat, vector<string> channels)
{
	cout<<"------------------------------------------------------------------------"<<endl;
	cout<<" Performing the fit using:                       "<< MCcat<< endl;
	cout<<" Considered fit categories:                      ";
	for(int i=0;i<channels.size();i++) cout<< channels.at(i) <<"  ";
	cout<<endl;
	cout<<" Is Barlow Beaston fit on:                       "<< BBeast<<endl;
	cout<<" Lc Ground state FF corrections activated:       "<< FFGstate<<endl;
	cout<<" Low Excited Lc state FF corrections activated:  "<< FFEstateL<<endl;
	cout<<" High Excited Lc state FF corrections activated: "<< FFEstateH<<endl;
	cout<<" Sweight corrections activated:                  "<< SWcorr << endl;
	cout<<"------------------------------------------------------------------------"<<endl;
}

map<string, Bool_t> DefineMapCategories(Bool_t FitIso, Bool_t FitKenr, Bool_t FitLcpi)
{
	map<string, Bool_t> FitMap;
	FitMap.insert(make_pair("Isolated", FitIso));
	FitMap.insert(make_pair("Kenriched",FitKenr));
	FitMap.insert(make_pair("Lcpipi",   FitLcpi));
	return FitMap;
}

vector<string> Channels2Fit(map<string,Bool_t> FitMap)
{
	vector<string> channels;
	map<string, Bool_t>::iterator itr;
	for (itr =  FitMap.begin(); itr !=  FitMap.end(); ++itr)
	{
		if(itr->second)
			channels.push_back(itr->first);
	}

	return channels;
}

void SetParametersIsolatedCh(SysFit &a, string ch)
{	
	a.SetStartParameter(ch,"Lb_Lcmunu",8.E5,3.E3,2.E6);
	a.SetStartParameter(ch,"Lb_LcDs-2body",14.E3,1E2,2.E5);
	a.SetStartParameter(ch,"Lb_LcDs-mbody",6.0E3,1E2,1.E5);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-2body",2.0E3,1E2,1.E5);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-2body",2.0E3,1E2,1.E5);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-mbody",2.0E3,1E2,1.E5);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-mbody",2.0E3,1E2,1.E5);
	a.SetStartParameter(ch,"Lb_Lc2625munu",1.8E5,1.E3,1.E6);
	a.SetStartParameter(ch,"Lb_Lc2593munu",0.8E5,1.E3,1.E6);
	a.SetStartParameter(ch,"MISID",3.1E4,10,1.E6);
	a.SetStartParameter(ch,"Combinatorial",4.0E4,10,1.E6);
	a.SetStartParameter(ch,"Lb_Lctaunu",0.05,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2625taunu",0.05,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2593taunu",0.05,0.,1);
}

void SetParametersKenrichedCh(SysFit &a, string ch)
{	
	a.SetStartParameter(ch,"Lb_Lcmunu",2E2,0,1.E3);
	a.SetStartParameter(ch,"Lb_LcDs-2body",3.E3,0,1.E4);
	a.SetStartParameter(ch,"Lb_LcDs-mbody",2.E3,500,2.E4);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-2body",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-2body",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-mbody",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-mbody",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2625munu",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2593munu",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"MISID",1.7E3,100,1.E4);
	a.SetStartParameter(ch,"Combinatorial",1.7E3,1.E2,5.E3);
	a.SetStartParameter(ch,"Lb_Lctaunu",0.,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2625taunu",0.,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2593taunu",0.,0.,1);
}

void SetParametersLcpipiCh(SysFit &a, string ch)
{	
	a.SetStartParameter(ch,"Lb_Lcmunu",2E2,0,1.E3);
	a.SetStartParameter(ch,"Lb_LcDs-2body",1.E3,0,1.E4);
	a.SetStartParameter(ch,"Lb_LcDs-mbody",1.E2,0,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-2body",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-2body",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2625Ds-mbody",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2593Ds-mbody",5.0E2,10,1.E4);
	a.SetStartParameter(ch,"Lb_Lc2625munu",5.0E3,10,2.E4);
	a.SetStartParameter(ch,"Lb_Lc2593munu",5.0E3,10,2.E4);
	a.SetStartParameter(ch,"MISID",5.E2,0,1.E3);
	a.SetStartParameter(ch,"Combinatorial",5.E1,0,5.E2);
	a.SetStartParameter(ch,"Lb_Lctaunu",0.,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2625taunu",0.05,0.,1);
	a.SetStartParameter(ch,"Lb_Lc2593taunu",0.05,0.,1);
}

string GetTemplateFileName(string ch, string MCcat)
{	
	string fname = "TemplateFiles/Histos_"+ch+"_"+MCcat+"_LbCorr";
	vector<string> suffix = {};
	if(FFGstate)
		suffix.push_back("_FFGstate");
	if(FFEstateL)
		suffix.push_back("_FFEstateL");
	if(FFEstateH)
		suffix.push_back("_FFEstateH");
	if(!FFGstate && !FFEstateL && !FFEstateH)
		suffix.push_back("_NoFFcorr");
	for (Int_t i=0; i<suffix.size();i++)
		fname+=suffix[i];
	fname += ".root";
	cout<<fname<<endl;
	
	return fname;
}

void RunFit(bool rebuild=true)
{
	// Check if the class is already loaded

	if (!TClass::GetDict("SysFit")||rebuild==true)
		gROOT->ProcessLine(".L SysFit.h++");


	string MCcat = "MCfull";
	if (MCTO) MCcat = "MCTrackerOnly";

	map<string, Bool_t> FitMap  = DefineMapCategories(FitIso, FitKenr, FitLcpi);
	vector<string> channels = Channels2Fit(FitMap);

	PrintFitConfiguration(MCcat, channels);

	//Define the fit
	SysFit a;

	//Set the fit configuration
	a.AllowBarlowBeaston(BBeast);
	a.CorrectFF_GS(FFGstate);
	a.CorrectFF_LES(FFEstateL);
	a.CorrectFF_HES(FFEstateH);
	a.CorrectSweights(SWcorr);
	a.SetMCcategory(MCcat);
	a.SelectChannel2fit(channels);

	//Set the start parameters of the fit
	for (int i=0;i<channels.size();i++)
	{
		string fname = GetTemplateFileName(channels[i], MCcat);
		a.SetTemplateFileName(channels[i],fname);
		if(channels[i]=="Isolated")
			SetParametersIsolatedCh(a, channels[i]);
		if(channels[i]=="Kenriched")
			SetParametersKenrichedCh(a, channels[i]);
		if(channels[i]=="Lcpipi")
			SetParametersLcpipiCh(a, channels[i]);
	}

	//Run the fit only
	RooStats::HistFactory::Measurement m = a.CreateMeasurement();
	RooWorkspace* w = a.CreateWorkspace(m);
    RooStats::ModelConfig* mc = a.CreateModel(w);
	RooFitResult *fitResult = a.Fit(mc, m, w, false,false);

	//Produce plot
    RooStats::HistFactory::Measurement m2 = a.CreateMeasurement();
    RooWorkspace* w2 = a.CreateWorkspace(m2);
    RooStats::ModelConfig* mc2 = a.CreateModel(w2);
    RooFitResult *fitResult2 = a.Fit(mc2, m2, w2, true, false);


}

