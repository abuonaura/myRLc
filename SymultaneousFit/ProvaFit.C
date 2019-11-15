#include "SysFit.h"
#include "SysFit.cpp"
#include <TROOT.h>
#include <vector>

//using namespace RLC;

void ProvaFit(Bool_t rebuild=true)
{
	cout<<rebuild<<endl;

	 // Check if the class is already loaded
	 if (!TClass::GetDict("SysFit")||rebuild==true)
	 {
    		gROOT->ProcessLine(".L SysFit.h++");
    		//gROOT->ProcessLine(".L SysFit.cpp++");
	 }
	 
	 TClass::GetClass("SysFit")->Print();

	SysFit a;


	//vector<string> ch2fit = {"Kenriched"};
	vector<string> ch2fit = {"Isolated"};
	a.SelectChannel2fit(ch2fit);
	vector<string> names = a.NameChannels();

	
	for (Int_t i =0; i<names.size();i++)
	{
		cout<<names[i]<<endl;
		string filename = string("RootFiles/Histos_")+names[i]+string(".root");
		cout<<filename<<endl;
		map<string,vector<Double_t>> startparams = a.GetStartParameters(names[i]);
		a.PrintStartParams(names[i],startparams);

		vector<string> category = a.GetCategory(startparams);
		for(Int_t j=0; j<category.size(); j++)
		{
			cout<<"category: "<<category[j]<<endl;
			if (category[j]=="pha")
				continue;
			cout<<a.GetHistoNormalisation(filename, category[j])<<endl;
			if(category[j]=="2charm-mbody")
			{
				a.ActivateShapeUncertainties(category[j],true);
				a.ActivateGaussConstraint(category[j],false);
			}
			else
			{
				a.ActivateShapeUncertainties(category[j],false);
				if(category[j]=="MISID"||category[j]=="Combinatorial")
					a.ActivateGaussConstraint(category[j],false);
				else
					a.ActivateGaussConstraint(category[j],false);
			}
			cout<<"--------------------------------------------"<<endl;
			cout<<"Is category : "<<category[j]<<"  shape uncertain? "<< a.IsShapeUncertain(category[j])<<endl;
			cout<<"Is category : "<<category[j]<<"  Gauss constrained? "<< a.IsGaussConstrained(category[j])<<endl;
			cout<<"--------------------------------------------"<<endl;
		}

		

	}
	
	RooStats::HistFactory::Measurement m = a.CreateMeasurement();
	RooWorkspace* w = a.CreateWorkspace(m);
	
	RooStats::ModelConfig* mc = a.CreateModel(w);
	RooFitResult *fit = a.Fit(mc, m, w);
	
	for (Int_t i =0; i<names.size();i++)
	{
		a.SetStartParameters(fit, names[i]);
		map<string,vector<Double_t>> newstartpars = a.GetStartParameters(names[i]);
		a.PrintStartParams(names[i],newstartpars);
	}

	//a.PrintStartParameters(
	//
/*
	a.AllowBarlowBeaston();	
	m = a.CreateMeasurement();
	w = a.CreateWorkspace(m);
	mc = a.CreateModel(w);
	RooFitResult* fit_new = a.Fit(mc, m, w);
*/
//	cout<<"---- "<<a.IsShapeUncertain("2charm-mbody")<<endl;
}
