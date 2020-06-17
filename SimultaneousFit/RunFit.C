#include "SysFit.h"
#include "SysFit.cpp"
#include <TROOT.h>
#include <vector>

void RunFit(Bool_t rebuild=true,string MCcat="MCfull",string FitType="Kenriched")
{
    cout<<rebuild<<endl;

     // Check if the class is already loaded
     if (!TClass::GetDict("SysFit")||rebuild==true)
     {
            gROOT->ProcessLine(".L SysFit.h++");
            //gROOT->ProcessLine(".L SysFit.cpp++");
     }

     TClass::GetClass("SysFit")->Print();

	//Set type of fit to be run
    SysFit a;
	a.SetMCcathegory(MCcat);
	if (FitType=="Kenriched")
		a.FitKenriched();
	if (FitType=="Isolated")
		a.FitIsolated();
	if (FitType=="Simultaneous")
		a.DoSimultaneousFit();

	//Get name of channels for the fit (Isolated/Kenriched)
	vector<string> names = a.NameChannels();

	for (Int_t i =0; i<names.size();i++)
    {	
		//File with templates
		string filename = string("RootFiles/Histos_")+names[i]+string("_")+MCcat+string(".root");
		//Get start parameters for fit values
		map<string,vector<Double_t>> startparams = a.GetStartParameters(names[i]);
		
		vector<string> category = a.GetCategory(startparams);
		for(Int_t j=0; j<category.size(); j++)
        {
            cout<<"category: "<<category[j]<<endl;
            if (category[j]=="pha")
                continue;
			
			//Do not use gauss constraints for any category
            a.ActivateGaussConstraint(category[j],false);
			
			//Set to false shape uncertainties (activate only for 2charm)
			a.ActivateShapeUncertainties(category[j],false);
			//Activate shape uncertainties for 2charm (groun/excited state)
			if(category[j]=="2charm-mbody")
            {
                a.ActivateShapeUncertainties(category[j],true);
            }
            else if(category[j]=="starDs-mbody")
            {
                a.ActivateShapeUncertainties(category[j],true);
			}
		}
	}

	//Perform fit
	RooStats::HistFactory::Measurement m = a.CreateMeasurement();
    RooWorkspace* w = a.CreateWorkspace(m);

    RooStats::ModelConfig* mc = a.CreateModel(w);
    RooFitResult *fit = a.Fit(mc, m, w);

}

