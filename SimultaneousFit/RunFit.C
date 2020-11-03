#include "SysFit.h"
#include "SysFit.cpp"
#include <TROOT.h>
#include <vector>

void RunFit(string MCcat,string FitType,Bool_t ffcorr=false, Bool_t rebuild=true)
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
	if (FitType=="Lcpipi")
		a.FitLcpipi();
	if (FitType=="Isolated")
		a.FitIsolated();
	if (FitType=="Simultaneous")
		a.DoSimultaneousFit();

	a.ActivateFFCorrections(ffcorr);

	a.Constrain2body(0.3);
	a.ConstrainMbody(0.5);

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
			else if(category[j]=="mu" || category[j]=="tau")
			{
				if (ffcorr)
					a.ActivateShapeUncertainties(category[j],true);
				else
					a.ActivateShapeUncertainties(category[j],false);
			}
		}
	}

	//Perform fit
	a.CorrectSweights(true);
	RooStats::HistFactory::Measurement m = a.CreateMeasurement();
    RooWorkspace* w = a.CreateWorkspace(m);

    RooStats::ModelConfig* mc = a.CreateModel(w);
    RooFitResult *fitResult = a.Fit(mc, m, w, false);

	//Produce plot
	a.CorrectSweights(false);
    RooStats::HistFactory::Measurement m1 = a.CreateMeasurement();
    RooWorkspace* w1 = a.CreateWorkspace(m1);

    RooStats::ModelConfig* mc1 = a.CreateModel(w1);
    RooFitResult *fitResult1 = a.Fit(mc1, m1, w1, true);

	for (Int_t i =0; i<names.size();i++)
    {
		string outfilename;
		if(FitType=="Simultaneous")
			outfilename = string("StoredTxtFitResults/StoredFitResults_")+names[i]+string("_")+MCcat+string("_Simultaneous.txt");
		else
			outfilename = string("StoredTxtFitResults/StoredFitResults_")+names[i]+string("_")+MCcat+string(".txt");
        a.StoreFitResults(outfilename,fitResult);
        a.CheckDiscrepancyWrtLastRLcValue(outfilename, FitType);
	}

}

