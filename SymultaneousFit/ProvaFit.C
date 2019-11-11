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
		}

	}

	a.Fit();


}
