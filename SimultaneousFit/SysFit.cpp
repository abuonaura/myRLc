#include "SysFit.h"
#include <TStyle.h>
#include <TRandom3.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace RooStats;
using namespace RooFit;
using namespace HistFactory;

SysFit::SysFit()
{
	alpha = 17;
	beta = 34;
	gamma = 98;
	alpha_s = 456;
	beta_s = 172;
	gamma_s = 23;
	BBeast = false;
	Twobodyconstraint=false;
	Mbodyconstraint=false;
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncmu_Isolated",{8.E5,3.E3,2.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nc2charm-2body_Isolated",{14.E3,1E2,2.E5}));
    start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nc2charm-mbody_Isolated",{6.0E3,1E2,1.E5}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcstarDs-2body_Isolated",{4.E3,1E2,1.E5}));
    start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcstarDs-mbody_Isolated",{4.E3,1E2,1.E5}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncstarmu_Isolated",{2.6E5,3.E4,1.E6}));
	//start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcLcpbar_Isolated",{2.6E5,3.E4,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcMISID_Isolated",{3.1E4,10,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcCombinatorial_Isolated",{4.0E4,1.E1,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nctau_Isolated",{0.05,0.,1.}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncstartau_Isolated",{0.05,0.,1.}));

	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncmu_Kenriched",{2E2,0,1.E3}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncstarmu_Kenriched",{1.E3,0,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncstartau_Kenriched",{0.,0.,1}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nctau_Kenriched",{0.,0.,1.}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nc2charm-2body_Kenriched",{3.E3,0,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nc2charm-mbody_Kenriched",{2.E3,500,2.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcstarDs-2body_Kenriched",{1.E3,10,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcstarDs-mbody_Kenriched",{1.E3,10,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcCombinatorial_Kenriched",{1.7E3,1.E2,5.E3}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcMISID_Kenriched",{1.7E3,100,1.E4}));

	start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Ncmu_Lcpipi",{2E2,0,1.E3}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Ncstarmu_Lcpipi",{1.E4,1.E3,2.E4}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Ncstartau_Lcpipi",{0.5,0.,1}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Nctau_Lcpipi",{0.,0,1.}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Nc2charm-2body_Lcpipi",{1.E3,0,1.E4}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("Nc2charm-mbody_Lcpipi",{1.E2,0,1.E4}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("NcstarDs-2body_Lcpipi",{1.E3,0,1.E4}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("NcstarDs-mbody_Lcpipi",{1.E3,0,1.E4}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("NcCombinatorial_Lcpipi",{5.E1,0,5.E2}));
    start_parameters["Lcpipi"].insert(pair<string, vector<Double_t>>("NcMISID_Lcpipi",{5E2,0,1.E3}));
}

void SysFit::FitIsolated()
{
	FitType = "Single";
	vector<string> ch2fit = {"Isolated"};
	SelectChannel2fit(ch2fit);
}

void SysFit::FitKenriched()
{
	FitType = "Single";
	vector<string> ch2fit = {"Kenriched"};
	SelectChannel2fit(ch2fit);
}

void SysFit::FitLcpipi()
{
    FitType = "Single";
    vector<string> ch2fit = {"Lcpipi"};
    SelectChannel2fit(ch2fit);
}

void SysFit::DoSimultaneousFit()
{
	FitType = "Simultaneous";
	vector<string> ch2fit = {"Isolated","Kenriched","Lcpipi"};
	//vector<string> ch2fit = {"Kenriched","Lcpipi"};
	SelectChannel2fit(ch2fit);
}

void SysFit::SetWeightGaussConstraint(string sample, Double_t w)
{
	weight[sample] = w;
}

Double_t SysFit::GetWeightGaussConstraint(string sample)
{
	Double_t w = weight[sample];
	return w;
}

void SysFit::ActivateShapeUncertainties(string sample, Bool_t value)
{
	ShapeUnc[sample] = value;
}

Bool_t SysFit::IsShapeUncertain(string sample)
{
	Bool_t v = ShapeUnc[sample];
	return v;
}

void SysFit::ActivateGaussConstraint(string sample, Bool_t value)
{
	GaussConstr[sample] = value;
}

Bool_t SysFit::IsGaussConstrained(string sample)
{
	Bool_t v = GaussConstr[sample];
	return v;
}

map<string,vector<Double_t>> SysFit::GetStartParameters(string ch_name)
{
	map<string,vector<Double_t>> params = start_parameters[ch_name];
	
	return params;
}



void SysFit::PrintStartParams(string channel, map<string,vector<Double_t>> start_parameters)
{
	cout<<"----  Paramameters channel:  "<<channel<<endl;
	// declaring iterators
	map<string, vector<Double_t>>::iterator it1 ;
	for (it1 = start_parameters.begin(); it1!=start_parameters.end(); it1++)
	{
		cout << it1->first << endl;
		cout<<start_parameters[it1->first][0]<<"  "<<start_parameters[it1->first][1]<<"   "<< start_parameters[it1->first][2]<<endl;

	}
}

vector<string> SysFit::GetParametersName(map<string,vector<Double_t>> parameters)
{       
	vector<string> param_name;
	// declaring iterators
	map<string, vector<Double_t>>::iterator it1 ;
	for (it1 =parameters.begin(); it1!=parameters.end(); it1++)
	{       
			param_name.push_back(string(it1->first));
	}

	return param_name;
}

vector<string> SysFit::GetCategory(map<string,vector<Double_t>> parameters)
{       
	vector<string> category;
	// declaring iterators
	map<string, vector<Double_t>>::iterator it1 ;
	for (it1 = parameters.begin(); it1!=parameters.end(); it1++)
	{      
		string name = string(it1->first);
		size_t f1 = name.find('_');
		category.push_back(name.substr(2, f1 - 2));
	}

	return category;
}

TString SysFit::GetComponentName(TString component)
{
	TString name = "";
	if (component.Contains("_mu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #mu #nu";
	else if(component.Contains("h_starmu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}* #mu #nu_{#mu}";
	else if(component.Contains("h_starDs-2body")) name = "#Lambda_{b} #rightarrow #Lambda_{c}* X_{c}";
	else if(component.Contains("h_starDs-mbody")) name = "#Lambda_{b} #rightarrow #Lambda_{c}* X_{c} X";
	else if(component.Contains("h_startau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}* #tau #nu_{#tau}";
	else if(component.Contains("h_Lcpbar_")) name = "#B #rightarrow #Lambda_{c} #bar{p} #mu_{#mu} #nu_{#tau}";
	else if(component.Contains("_tau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #tau #nu_{#tau}";
	else if(component.Contains("h_2charm-2body")) name = "#Lambda_{b} #rightarrow #Lambda_{c} X_{c}";
	else if(component.Contains("h_2charm-mbody")) name = "#Lambda_{b} #rightarrow #Lambda_{c} X_{c} X";
	else if(component.Contains("h_MISID")) name = "MISID";
	else if(component.Contains("h_Combinatorial")) name =  "Combinatorial";

	return name;
}

Int_t SysFit::GetComponentColor(TString component)
{       
	Int_t color=-99;
	if (component.Contains("_mu_")) color = kBlue;
	else if(component.Contains("h_starmu_")) color = kViolet;
	else if(component.Contains("h_startau_")) color = kMagenta-9;
	else if(component.Contains("h_Lcpbar_")) color = kViolet+7;
	else if(component.Contains("h_starDs-2body")) color = kAzure+2;
	else if(component.Contains("h_starDs-mbody")) color = kAzure+10;
	else if(component.Contains("_tau_"))color = kRed;
	else if(component.Contains("h_2charm-2body")) color = kGreen;
	else if(component.Contains("h_2charm-mbody")) color = kGreen+3;
	else if(component.Contains("h_MISID")) color = kYellow;
	else if(component.Contains("h_Combinatorial")) color = kOrange+8;
	return color;
}

TString SysFit::GetFitVarName(TString title)
{
	TString name = "";
	if (title.Contains("q")) name = "q2";
	else if (title.Contains("E")) name = "El";
	else if (title.Contains("M")) name = "Mmiss";
	return name;
}

Double_t SysFit::GetHistoNormalisation(string inputFile, string type)
{
	cout<<" InputFile name: " << inputFile<<endl;
	cout<<" type: " <<type<<endl;
	//Read root file with histograms of samples
	TFile histFile(inputFile.c_str());
	TH1 *h_temp;
	Bool_t ffcorr = GetFFcorrectionValue();
	//Get the histo from the file and save in htemp
	
	if (type=="mu" || type=="tau")
    {
		if(ffcorr)
			histFile.GetObject(("h_w_"+type+"_mean").c_str(),h_temp);
		else
			histFile.GetObject(("h_"+type).c_str(),h_temp);
	}
	else
		histFile.GetObject(("h_"+type).c_str(),h_temp);
	//Verify the histo is non null
	assert(h_temp!=NULL);

	//Get Normalisation factor to 1
	Double_t mcNorm = 1./h_temp->Integral();

	cout<<" Integral_"+type+" = "<<h_temp->Integral()<<endl;
	cout<<" Normalisation factor "+type+" = " << mcNorm <<endl;

	return mcNorm;
}


void SysFit::CreateSweightCorrectHistos(string inputFile, string inputFile2,string newfilename, vector<string> category)
{
    //Read root file with histograms of samples
    Bool_t ffcorr = GetFFcorrectionValue();
    TFile *oldhistFile = new TFile(inputFile.c_str(),"READ");
    TFile *oldhistFile2 = new TFile(inputFile2.c_str(),"READ");

    TH3* h_data;
    TH3* sweightcorrhisto;

    h_data = (TH3*)oldhistFile->Get("h_data");

    TFile *newhistFile = new TFile(newfilename.c_str(),"Recreate");
    cout<<"new file name: "<<newfilename<<endl;
    sweightcorrhisto = (TH3*)h_data->Clone("h_data_swcorr");
    sweightcorrhisto->SetDirectory(0);
    sweightcorrhisto->Reset();

    double corr_num, corr_den;

    for (unsigned int xbin_idx=1; xbin_idx<=h_data->GetXaxis()->GetNbins(); xbin_idx++)
    {
        for (unsigned int ybin_idx=1; ybin_idx<=h_data->GetYaxis()->GetNbins(); ybin_idx++)
        {
            for (unsigned int zbin_idx=1; zbin_idx<=h_data->GetZaxis()->GetNbins(); zbin_idx++)
            {
                corr_num = h_data->GetBinContent(xbin_idx,ybin_idx,zbin_idx);
                corr_den = h_data->GetBinError(xbin_idx,ybin_idx,zbin_idx)*h_data->GetBinError(xbin_idx,ybin_idx,zbin_idx);
                if(corr_den==0.)
                {
                    corr_num = 1.;
                    corr_den = 1.;
                }
                //cout<<corr_num<<"   "<<corr_den<<"    "<<corr_num/corr_den<<endl;
                 sweightcorrhisto->SetBinContent(xbin_idx,ybin_idx,zbin_idx,corr_num/corr_den);
                 sweightcorrhisto->SetBinError(xbin_idx,ybin_idx,zbin_idx,0.);
            }
        }
    }
    sweightcorrhisto->Write();
    h_data->Sumw2();
    h_data->Multiply(sweightcorrhisto);
    h_data->Write();

	for(Int_t j=0; j<category.size(); j++)
    {
        //Get the histo from the file and save in htemp
        if (category[j]!="2charm-mbody"&&category[j]!="starDs-mbody"&&category[j]!="mu"&&category[j]!="tau")
        {
            TH3* h_temp;
            oldhistFile->GetObject(("h_"+category[j]).c_str(),h_temp);
            if(h_temp==NULL)
                cout<<"-------- hetemp Null! ----"<<endl;

            h_temp->Sumw2();
            h_temp->Multiply(sweightcorrhisto);
            h_temp->Write();
        }
        else if(category[j]=="2charm-mbody" || category[j]=="starDs-mbody")
        {
            std::vector<TH3*>* h_temp = new std::vector<TH3*>{};
            vector<string> suffix = {"1pl","1ml","1pq","1mq"};
            h_temp->push_back((TH3*)oldhistFile->Get(("h_"+category[j]).c_str()));
            for(Int_t a = 0; a<suffix.size();a++)
            {
                h_temp->push_back((TH3*)oldhistFile->Get(("h_"+category[j]+"_"+suffix[a]).c_str()));
            }
            for(Int_t a = 0; a<h_temp->size();a++)
            {
                cout<<a<<endl;
                if(h_temp->at(a)==NULL)
                    cout<<"-------- hetemp Null! ----"<<endl;
                h_temp->at(a)->Sumw2();
                h_temp->at(a)->Multiply(sweightcorrhisto);
                h_temp->at(a)->Write();
            }

        }
		else if(category[j]=="mu" || category[j]=="tau")
        {
            if(ffcorr)
            {
                std::vector<TH3*>* h_temp = new std::vector<TH3*>{};
                vector<string> suffix = {"mean","max","min"};
                for(Int_t a = 0; a<suffix.size();a++)
                {
                    h_temp->push_back((TH3*)oldhistFile2->Get(("h_w_"+category[j]+"_"+suffix[a]).c_str()));
                    if(h_temp->at(a)==NULL)
                        cout<<"-------- hetemp Null! ----"<<category[j]<<"  "<<suffix[a]<<endl;
                    h_temp->at(a)->Sumw2();
                    h_temp->at(a)->Multiply(sweightcorrhisto);
                    h_temp->at(a)->Write();
                }
            }
            else
            {
                TH3* h_temp;
                oldhistFile->GetObject(("h_"+category[j]).c_str(),h_temp);
                h_temp->Sumw2();
                h_temp->Multiply(sweightcorrhisto);
                h_temp->Write();
            }


        }
    }
    newhistFile->Close();
    oldhistFile->Close();
    oldhistFile2->Close();

}

void SysFit::AddSample(string type, string inputFile, bool shapeUncert, bool GaussConstraints, const bool BBeast, Channel** chan, vector<Double_t> params, vector<Double_t> mu_params={0})
{
	//Define the sample
    Sample sample;

    //set the file from where to take the MC normalisation
    string MCnormFile = inputFile;

    //set the file from where to take the template
	string TemplateFile = inputFile;
    if(CorrectSweightsValue()==true)
      TemplateFile = GetSweightHistFileName();
    sample.SetInputFile(TemplateFile);

    cout<<"---------------------------"<<endl;
    cout<<"File from which the normalisation is read: "<<MCnormFile<<endl;
    cout<<"File from which the template is read: "<<TemplateFile<<endl;
    cout<<"---------------------------"<<endl;

    //get the normalisation factor to 1
    Double_t mcNorms = GetHistoNormalisation(MCnormFile, type);

	//Get the name of the channel
    string chname =(*chan)->GetName();
    string category;
    if (chname.find("Isolated") != std::string::npos) category = "Isolated";
    if (chname.find("Kenriched") != std::string::npos) category = "Kenriched";
    if (chname.find("Lcpipi") != std::string::npos) category = "Lcpipi";

	if(!shapeUncert)
	{
			sample.SetName("h_"+type);
            sample.SetHistoName("h_"+type);
	}
	else
	{
		if(type=="2charm-mbody"||type=="starDs-mbody")
		{
			sample.SetName("h_"+type);
			sample.SetHistoName("h_"+type);
			sample.AddHistoSys(type+"_variation","h_"+type+"_1ml",TemplateFile,"", "h_"+type+"_1pl",TemplateFile,"");
			sample.AddHistoSys(type+"_quadratic_variation","h_"+type+"_1mq",TemplateFile,"", "h_"+type+"_1pq",TemplateFile,"");
		}
		if(type=="mu"||type=="tau")
		{
			sample.SetName("h_w_"+type+"_mean");
			sample.SetHistoName("h_w_"+type+"_mean");
			sample.AddHistoSys(type+"_shape_unc","h_w_"+type+"_min",TemplateFile,"","h_w_"+type+"_max",TemplateFile,"");
		}
	}

	//set a sample to be "normalised by theory" (its normalisation scales with luminosity)
    sample.SetNormalizeByTheory(kFALSE);

	//Normalise sample to 1 (FOR ALL samples)
	sample.AddNormFactor("mcNc"+type+"_"+category,mcNorms,1e-20,1.);

	//Add additional normalisation factors
	if(type=="tau")
	{
		sample.AddNormFactor("Ncmu_"+category,mu_params[0],mu_params[1],mu_params[2]); //normalisation to the start params of the Lc muon sample
		sample.AddNormFactor("RLc"+type+"_"+category, params[0],params[1],params[2]); //times the value of RLc
	}
	else if(type=="startau")
    {
        sample.AddNormFactor("Ncstarmu_"+category,mu_params[0],mu_params[1],mu_params[2]); //normalisation to the star params of the Lambda_c* muon sample
        sample.AddNormFactor("RLc"+type+"_"+category, params[0],params[1],params[2]); //times the value of RLc star
    }
	else if(type=="starDs-2body"&&Twobodyconstraint==true)
	{
		sample.AddNormFactor("Nc2charm-2body_"+category,params[0],params[1],params[2]); //normalisation to the star params of the Lambda_cXc-2body sample
		sample.AddNormFactor("NcRatio2body",ratio2body,1e-9,1.);
	}
	else if(type=="starDs-mbody"&&Mbodyconstraint==true)
	{
		sample.AddNormFactor("Nc2charm-mbody_"+category,params[0],params[1],params[2]); //normalisation to the star params of the Lambda_cXc-2body sample
		sample.AddNormFactor("NcRatioMbody",ratioMbody,1e-9,1.);
	}
	else
		sample.AddNormFactor("Nc"+type+"_"+category, params[0],params[1],params[2]); //normalisation to starting parameters (for all other samples)

	//----- Add Gaussian Constraint for MISID/Combinatorial for control fit part
    if (GaussConstraints && chname=="RLc_kinematic_Kenriched")
    {
        Double_t w = GetWeightGaussConstraint(type);
        sample.AddOverallSys("NormUnc_"+type+"_"+category,1.,w);
    }

    if(BBeast) sample.ActivateStatError();
    (*chan)->AddSample(sample);
}

RooStats::ModelConfig* SysFit::SetChannelConstants(RooStats::ModelConfig *mc, string channel)
{
	if (channel=="Kenriched")
	{
		((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel).c_str())))->setVal(1E-15);
		((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLcstartau_"+channel).c_str())))->setVal(1E-15);
		((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLcstartau_"+channel).c_str())))->setConstant(kTRUE);
		if (IsGaussConstrained("MISID"))
			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcMISID_"+channel).c_str())))->setConstant(kTRUE);
		if (IsGaussConstrained("Combinatorial"))
			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcCombinatorial_"+channel).c_str())))->setConstant(kTRUE);
	}
	if (channel=="Lcpipi")
    {   
        ((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel).c_str())))->setVal(1E-15);
        ((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel).c_str())))->setConstant(kTRUE);
    }
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcmu_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNctau_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarmu_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstartau_"+channel).c_str())))->setConstant(kTRUE);
	if((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcLcpbar_"+channel).c_str())))
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcLcpbar_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarDs-mbody_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarDs-2body_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm-mbody_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm-2body_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcMISID_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcCombinatorial_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find("NcRatio2body")))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find("NcRatioMbody")))->setConstant(kTRUE);
	
	return mc;
}

RooStats::ModelConfig* SysFit::FixYields(RooStats::ModelConfig *mc, string sample, string channel)
{
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("Nc"+sample+"_"+channel).c_str())))->setConstant(kTRUE);
	return mc;
}

RooStats::ModelConfig* SysFit::SetAlphaStartingPoints(RooStats::ModelConfig *mc)
{
	((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_2charm-mbody_variation")))->setVal(0.5);
	((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_starDs-mbody_variation")))->setVal(-0.2);
	((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_2charm-mbody_quadratic_variation")))->setVal(0.4);
	((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_starDs-mbody_quadratic_variation")))->setVal(-0.3);
	return mc;
}

RooStats::HistFactory::Measurement SysFit::CreateMeasurement()
{       

	//Define the measurement
	Measurement meas("RLc","Rb");
	//Define where to save the workspace
	meas.SetOutputFilePrefix("results/RLc");
	//Tell histfactory not to run the fit itself
	meas.SetExportOnly(kTRUE);
	//Set the luminosity
	meas.SetLumi(1.0);
	meas.SetLumiRelErr(0.05);

	//Define a channel for each category
	vector <string> channel_names = NameChannels();	
	TString MCcat = GetMCcathegory();
	Int_t nchannels = channel_names.size();
	vector <Channel*> channels;

	//Understand if FF corrections are or not activated
	Bool_t ffcorr = GetFFcorrectionValue();
	cout<<endl;
	cout<<" FF corrections: "<<ffcorr<<endl;
	cout<<endl;
	
	for(Int_t i=0; i<nchannels; i++)
	{
		channels.push_back(new Channel((string("RLc_kinematic_")+channel_names[i]).c_str()));
		channels[i]->SetStatErrorConfig(1e-5, "Poisson");
		map<string,vector<Double_t>> start_param = GetStartParameters(channel_names[i]);
		vector<string> category = GetCategory(start_param);
		vector<string> param_names = GetParametersName(start_param);
		string filename = string("RootFiles/Histos_")+channel_names[i]+string("_")+string(MCcat)+string(".root");
		string filename1 = string("RootFiles/DemoHistosLattice_")+channel_names[i]+string("_")+string(MCcat)+string(".root");

		if(CorrectSweightsValue())
        {
            //Creating name for histograms with sweight correction
            string start = filename;
            string from = ".root";
            string to = "_swcorr.root";
            size_t start_pos = start.find(from);
            string newfilename = start.replace(start_pos, from.length(), to);
            SetSweightHistFileName(newfilename);

            CreateSweightCorrectHistos(filename, filename1, newfilename,category);
        }



		//Add the samples to the channels
		for(Int_t j=0; j<category.size(); j++)
		{
			cout<<category[j]<<endl;
			cout<<IsShapeUncertain(category[j])<<endl;
			if(category[j]=="tau")//NB:: To modify the 0 adding ShapeUnc
			{
				if(ffcorr)
					AddSample(category[j],filename1,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]],start_param[string("Ncmu_")+channel_names[i]]);
				else
					AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]],start_param[string("Ncmu_")+channel_names[i]]);
			}
			else if(category[j]=="mu")//NB:: To modify the 0 adding ShapeUnc
			{
				if(ffcorr)
					AddSample(category[j],filename1,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]]);
				else
					AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]]);
			}
			else if(category[j]=="startau")
				AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]],start_param[string("Ncstarmu_")+channel_names[i]]);
			else if(category[j]=="starDs-2body"&& Twobodyconstraint==true)
				AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[string("Nc2charm-2body_")+channel_names[i]]);
			else if(category[j]=="starDs-mbody"&& Mbodyconstraint==true)
				AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[string("Nc2charm-mbody_")+channel_names[i]]);
			else
				AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]]);
		}

		channels[i]->SetData("h_data",filename);
		Data ch_data = channels[i]->GetData();
		ch_data.SetName("h_data_"+channel_names[i]);
		cout<<")****** "<<ch_data.GetName()<<endl;
		meas.AddChannel(*channels[i]);
		//Define the parameter of interest
		meas.SetPOI("RLctau_"+channel_names[i]);
		meas.SetPOI("RLcstartau_"+channel_names[i]);
 
	}
	cout << endl;
	cout << "----------------------   Collecting histograms   ------------------------------" << endl;
	meas.CollectHistograms();
	cout << "-------------------------------------------------------------------------------" << endl;
	// Print to the screen a text representation of the model
	// just for minor debugging
	cout << "-------------------------   Printing Tree   ------------------------------" << endl;
	meas.PrintTree();
	cout << "-------------------------------------------------------------------------------" << endl;

	return meas;
}


RooWorkspace* SysFit::CreateWorkspace(RooStats::HistFactory::Measurement meas)
{
	//build a workspace with a pdf and a modelconfig
	RooWorkspace *w = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
	cout<<"--------------------- Workspace created -------------------"<<endl;
	//w->Print();
	//cout<<"-----------------------------------------------------------"<<endl;
	return w;
}


RooStats::ModelConfig* SysFit::CreateModel(RooWorkspace* w)
{
	//Get the model manually
	RooStats::ModelConfig *mc = (RooStats::ModelConfig*) w->obj("ModelConfig");
	return mc;
}

RooFitResult* SysFit::Fit(RooStats::ModelConfig* mc, RooStats::HistFactory::Measurement meas, RooWorkspace *w, Bool_t plotonly, Bool_t refit)
{
	//Understand if FF corrections are or not activated
	Bool_t ffcorr = GetFFcorrectionValue();	
	vector <string> channel_names = NameChannels();
	Int_t nchannels = channel_names.size();

	RooSimultaneous *model_ = (RooSimultaneous*)mc->GetPdf();

	//Fix the MC normalisation of the histos
	cout << "-------------------------------   Fixing parameters   ------------------------------" << endl;
	//mc->GetNuisanceParameters()->Print();
	((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);
	for(Int_t i=0; i<nchannels; i++)
		SetChannelConstants(mc,channel_names[i]);

	cout << "------------------------------------------------------------------------------------" << endl;


	cout << "------------------------------Setting Alpha Starting Values-------------------------" << endl;
	SetAlphaStartingPoints(mc);
	
	cout << "------------------Printig fit nuisance parameters starting values ------------------" << endl;
	//Print the values of the nuisance parameters befor the fit
	RooRealVar *tempVar; //dummy pointer
    TIterator *iterpars = mc->GetNuisanceParameters()->createIterator();
    while((tempVar = (RooRealVar*) iterpars->Next())) {
      TString theName = tempVar->GetName();
	  Double_t value = tempVar->getValV();
	  cout<<" Parameter "<<theName<<" starting from "<<value<<endl;
	}
	cout << "------------------------------------------------------------------------------------" << endl;

	
	if(plotonly||refit)
	{
	for(Int_t i=0; i<nchannels;i++)
        {
            string resultfname = "FitResultsUnblinded_"+channel_names[i]+".txt";
            LoadFitResults(resultfname, mc, channel_names[i]);
			cout << "------------------Printig fit nuisance parameters values (loaded) ------------------" << endl;
			//Print the values of the nuisance parameters befor the fit
			RooRealVar *tempVar; //dummy pointer
			TIterator *iterpars = mc->GetNuisanceParameters()->createIterator();
			while((tempVar = (RooRealVar*) iterpars->Next())) 
			{
				TString theName = tempVar->GetName();
				Double_t value = tempVar->getValV();
				cout<<" Parameter "<<theName<<" set to "<<value<<endl;
			}
			cout << "------------------------------------------------------------------------------------" << endl;
		}
	}


	//Tell the name of the obervables
	RooArgSet *obs = (RooArgSet*) mc->GetObservables();
	//obs->Print();

	std::vector<RooRealVar*> x_vector;
	std::vector<RooRealVar*> y_vector;
	std::vector<RooRealVar*> z_vector;

	for(Int_t i=0; i<nchannels;i++)
	{
		x_vector.push_back((RooRealVar*) obs->find(("obs_z_RLc_kinematic_"+channel_names[i]).c_str()));
		x_vector[i]->setUnit("[GeV^{2}/c^{4}]");
		y_vector.push_back((RooRealVar*) obs->find(("obs_y_RLc_kinematic_"+channel_names[i]).c_str()));
		y_vector[i]->setUnit("[MeV/c^{2}]");
		z_vector.push_back((RooRealVar*) obs->find(("obs_x_RLc_kinematic_"+channel_names[i]).c_str()));
		z_vector[i]->setUnit("[GeV^{2}/c^{4}]");
	}



	int num_channels = meas.GetChannels().size();
	//cout<< " num_channels = "<<num_channels<<endl;
	//cout<<endl;

	RooCategory *idx = (RooCategory*) obs->find("channelCat");
	RooAbsData *data = (RooAbsData*) w->data("obsData");

	int  nPlotsMax = 1000;
	RooDataSet* simData=NULL;
	RooSimultaneous* simPdf = NULL;
	if(strcmp(mc->GetPdf()->ClassName(),"RooSimultaneous")==0){
		cout <<"Is a simultaneous PDF"<<endl;
		simPdf = (RooSimultaneous *)(mc->GetPdf());
	} else {
		cout <<"Is not a simultaneous PDF"<<endl;
	}

	RooAbsPdf *pdf = model_->getPdf(idx->getLabel());

	RooSimultaneous *model = new RooSimultaneous("simPdf_modified","simPdf_modified",*idx);
	RooAbsPdf* pdf_;
	RooCustomizer* cust;

	//Put together the alphas from the FF corrections
	RooRealVar *alpha_mu_shape_unc;
	RooRealVar *alpha_tau_shape_unc;
	RooRealVar* alpha = new RooRealVar("alpha","alpha", 0, -3., 3.);
    if (ffcorr)
    {
        alpha_mu_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_mu_shape_unc");
        alpha_tau_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_tau_shape_unc");
    }

	//Put together the alphas from the mbody shape variations for LcXc and Lc*Xc
	RooRealVar *alphal_2charm;
    RooRealVar *alphal_starDs;
    RooRealVar *alphaq_2charm;
    RooRealVar *alphaq_starDs;
    RooRealVar* alpha_linear = new RooRealVar("alpha_linear","alpha_linear", 0, -3., 3.);
    RooRealVar* alpha_quadratic = new RooRealVar("alpha_quadratic","alpha_quadratic", 0, -3., 3.);
    alphal_2charm = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_2charm-mbody_variation");
    alphal_starDs = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_starDs-mbody_variation");
    alphaq_2charm = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_2charm-mbody_quadratic_variation");
    alphaq_starDs = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_starDs-mbody_quadratic_variation");

	cout<<"-----FF corr: "<<ffcorr<<endl;

	for (Int_t i =0;i< nchannels;i++)
	{
		//cout<<i<<endl;
		idx->setIndex(i);
		pdf_=model_->getPdf(idx->getLabel());
		cust = new RooCustomizer(*pdf_,"cust");
		/*
         * ---------- Make FF alpha parameters between tau and mu the same but different for the channels
        RooRealVar* alpha = new RooRealVar((string("alpha_")+channel_names[i]).c_str(),(string("alpha_")+channel_names[i]).c_str(), 0, -100., 100.);
        if (ffcorr)
        {
            alpha_mu_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find((string("alpha_mu_shape_unc_RLc_kinematic_")+channel_names[i]).c_str());
            alpha_tau_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find((string("alpha_tau_shape_unc_RLc_kinematic_")+channel_names[i]).c_str());
        }
        */

	if(plotonly==false)
	{


		if(ffcorr)
		{
			cust->replaceArg(*alpha_mu_shape_unc,*alpha);
			cust->replaceArg(*alpha_tau_shape_unc,*alpha);
		}
		//Set the linear alphas for the LcXc and LcXcStar templates to be one parameter
        cust->replaceArg(*alphal_2charm,*alpha_linear);
        cust->replaceArg(*alphal_starDs,*alpha_linear);
        //Set the quadratic alphas for the LcXc and LcXcStar templates to be one parameter
        cust->replaceArg(*alphaq_2charm,*alpha_quadratic);
        cust->replaceArg(*alphaq_starDs,*alpha_quadratic);
	}

		model->addPdf((RooAbsPdf&)*cust->build(),idx->getLabel());
	}

	if(plotonly==false)
	{
		cout<<" ------------------------- Starting fitting procedure ------------------------"<<endl;
		RooStats::HistFactory::HistFactorySimultaneous* model_hf = new RooStats::HistFactory::HistFactorySimultaneous(*model);

		RooAbsReal* nll_hf;
        if(refit==false) nll_hf	= model_hf->createNLL(*data);
		if(refit==true) nll_hf = model_hf->createNLL(*data,RooFit::Offset(kTRUE));//,RooFit::NumCPU(8));

		RooMinuit* minuit_hf = new RooMinuit(*nll_hf);
		RooArgSet *temp = new RooArgSet();

		minuit_hf->setErrorLevel(0.5);
		minuit_hf->setStrategy(2);
		minuit_hf->fit("smh");
		
		RooPlot* frame = y_vector[0]->frame();
		nll_hf->plotOn(frame) ;
		TCanvas *c =new TCanvas();
		frame->Draw();

		RooFitResult *fitResult=minuit_hf->save("TempResult","TempResult");

		for(Int_t i=0; i<nchannels;i++)
		{
			string resultfname = "FitResultsUnblinded_"+channel_names[i]+".txt";
			SaveFitResults(resultfname, fitResult);
		}


		for(Int_t i=0; i<nchannels;i++)
		{
			//if(channel_names[i]=="Isolated")
			blindResult(fitResult,channel_names[i]);
		}

		//cout <<"-------------------------- PRINT VERBOSE RESULTS -----------------------------" <<endl;
		//Verbose printing: Basic info, values of constant parameters, initial and
		// final values of floating parameters, global correlations
		//fitResult->Print("V");
		//cout <<"------------------------------------------------------------------------------" <<endl;

		//Summary printing: Basic info plus final values of floating fit parameters
		cout <<"-------------------------- PRINT SUMMARY RESULTS ------------------------------" <<endl;
		fitResult->Print();
		cout <<"------------------------------------------------------------------------------" <<endl;
	
		/*
		for(Int_t i=0; i<nchannels;i++)
    {
        idx->setIndex(i);
        PlotFitVariablesInFit(x_vector[i],"M_{miss}^{2}",y_vector[i],"E_{l}",z_vector[i],"q^{2}",data,model,idx,channel_names[i]);
    }
	*/
		return fitResult;
	}
	else
	{
		// Creation of the data datahists for plotting, and the inverse of the sWeight correction for the model.
		std::vector<TH3D> plotting_TH3D_vector;
		std::vector<RooDataHist> plotting_datahist_vector;

		for(Int_t i=0; i<nchannels;i++)
		{
			string filename = string("RootFiles/Histos_")+channel_names[i]+string("_")+string(GetMCcathegory())+string(".root");

			TFile* datatemplatefile = new TFile((filename).c_str(),"read");

			plotting_TH3D_vector.push_back(*((TH3D*) (datatemplatefile->Get("h_data"))->Clone(("TH1"+channel_names[i]).c_str())));
			datatemplatefile->Close();
			plotting_datahist_vector.push_back(RooDataHist(("plotting_datahist_channel"+channel_names[i]).c_str(),("plotting_datahist_channel"+channel_names[i]).c_str(),RooArgList(*z_vector[i],*y_vector[i],*x_vector[i]),&plotting_TH3D_vector[i]));
		}

		std::map<std::string,RooDataHist*> datahist_map_plotting;
		for (unsigned int i=0; i<nchannels; i++)
		{
			idx->setIndex(i);
			datahist_map_plotting.insert(std::pair<std::string,RooDataHist*>(std::string(idx->getLabel()),&plotting_datahist_vector[i]));
		}

		for(Int_t i=0; i<nchannels;i++)
		{ 
			idx->setIndex(i);
			PlotFitVariables(x_vector[i],"M_{miss}^{2}",y_vector[i],"E_{l}",z_vector[i],"q^{2}",datahist_map_plotting,model,idx,channel_names[i]);
		}

		return 0;

	}


	/*
	for(Int_t i=0; i<nchannels;i++)
    { 
        idx->setIndex(i);
		PlotFitVariables(x_vector[i],"M_{miss}^{2}",y_vector[i],"E_{l}",z_vector[i],"q^{2}",data,model,idx,channel_names[i]);
	}

	// Access list of final fit parameter values
	if(fitResult->status()==0)
		cout<<" Fit converged "<<endl;
	else
		cout<<"!!!! ATTENTION !!!!    Fit NON CONVERGING!  "<<endl;

	cout << "final value of floating parameters" << endl ;
	fitResult->floatParsFinal().Print("s") ;
	*/


	
	//RooRealVar* RLc_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLtau");
	//Double_t fitOut[] = {RLc_fitresult->getVal(),RLc_fitresult->errorVar()->getVal()};

}


void SysFit::blindResult(RooFitResult *fitResult,string name_suffix)
{
    RooRealVar* RLc_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLctau_"+TString(name_suffix));
	if(RLc_fitresult)
	{
		RLc_fitresult->setRange(-999999,999999);
		Double_t true_val = RLc_fitresult->getVal();
		Double_t error =  RLc_fitresult->getError();
		cout<<endl;
		cout<<endl;
		cout<<"------------------------------------------"<<endl;
		if (true_val>0 )
		{
			cout<< "DON'T WORRY! RLc is positive "<<endl;
			if ((true_val - 2*error)<0)
				cout<<" the value is however consistent with 0" <<endl;
			else
				cout<<" and the value is not consistent with 0" <<endl;
		}
		if (true_val<=0 )
		{
			cout<< "WORRY! RLc is either negative or zero! "<<endl;
			if ((true_val + 2*error)> 0)
				cout<<" and the value is consistent with 0" <<endl;
			else
				cout<<" and the value is not consistent with 0" <<endl;
			cout<<"Unblinded result: "<< true_val<<" +/- "<<  RLc_fitresult->getError()<<endl;
		}

		cout<<"------------------------------------------"<<endl;
		cout<<endl;
		cout<<endl;
		RLc_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha).Uniform(0,100)))*TRandom3(beta).Uniform(0,100)*true_val+TRandom3(gamma).Uniform(0,100));
	}
	 RooRealVar* RLcStar_fitresult;
   	
	 RLcStar_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLcstartau_"+TString(name_suffix));
	if(RLcStar_fitresult)
	{
		RLcStar_fitresult->setRange(-999999,999999);
		Double_t true_val_star=RLcStar_fitresult->getVal();
		Double_t error_star =  RLcStar_fitresult->getError();
		cout<<"------------------------------------------"<<endl;
		if (true_val_star>0 )
		{
			cout<< "DON'T WORRY! RLc* is positive "<<endl;
			if ((true_val_star - 2*error_star)<0)
                cout<<" the value is however consistent with 0" <<endl;
            else
                cout<<" and the value is not consistent with 0" <<endl;
			//cout<<"Unblinded result: "<< true_val_star<<" +/- "<<  RLcStar_fitresult->getError()<<endl;
		}
		if (true_val_star<=0 )
		{
			cout<< "WORRY! RLc* is either negative or zero! "<<endl;
			if ((true_val_star + 2*error_star)> 0)
                cout<<" and the value is consistent with 0" <<endl;
            else
                cout<<" and the value is not consistent with 0" <<endl;
			cout<<"Unblinded result: "<< true_val_star<<" +/- "<<  RLcStar_fitresult->getError()<<endl;
		}
		cout<<"------------------------------------------"<<endl;
		cout<<endl;
		cout<<endl;
		RLcStar_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha_s).Uniform(0,100)))*TRandom3(beta_s).Uniform(0,100)*RLcStar_fitresult->getVal()+TRandom3(gamma_s).Uniform(0,100));
	}
}

RooStats::ModelConfig* SysFit::LoadFitResults(string fname, RooStats::ModelConfig* mc, string channelName)
{
	ifstream infile;
	infile.open(fname);
	string search="RLctau_"+channelName;
	cout<< "string search: "<<search<<endl;
	string search2="RLcstartau_"+channelName;
	cout<< "string search: "<<search2<<endl;
	string search3 = "alpha_linear";
	string search4 = "alpha_quadratic";
	string search5 = "alpha";
	string name;
	Double_t value=0., error=0.;
	while(infile>>name>>value>>error)
	{
		if(name.compare(search)!=0 && name.compare(search2)!=0 && name.compare(search3)!=0 && name.compare(search4)!=0&& name.compare(search5)!=0)
		{	
			cout<<name<<endl;
			if((RooRealVar*)(mc->GetNuisanceParameters()->find(name.c_str())))
				((RooRealVar*)(mc->GetNuisanceParameters()->find(name.c_str())))->setVal(value);
			else
			{
				cout<<" Looking at channel: "<< channelName<<" parameter "<< name<< " not found. Leaving as it is"<<endl;
			}
		}
		else if (name.compare(search)==0 || name.compare(search2)==0)
			((RooRealVar*)(mc->GetParametersOfInterest()->find(name.c_str())))->setVal(value);
		else if(name.compare(search3)==0)
		{	
			((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_2charm-mbody_variation")))->setVal(value);
			((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_starDs-mbody_variation")))->setVal(value);
		}
		else if(name.compare(search4)==0)
		{	
			((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_2charm-mbody_quadratic_variation")))->setVal(value);
			((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_starDs-mbody_quadratic_variation")))->setVal(value);
		}
		else if(name.compare(search5)==0)
		{
			cout<< "   "<<name<<endl;
			((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_mu_shape_unc")))->setVal(value);
            ((RooRealVar*)(mc->GetNuisanceParameters()->find("alpha_tau_shape_unc")))->setVal(value);
		}
	}
	return mc;
}

void SysFit::SaveFitResults(string fname, RooFitResult *fitResult)
{
	ofstream outfile;
	outfile.open(fname,std::ios::trunc);
	auto pars = fitResult->floatParsFinal();
	RooRealVar *p;
	for (Int_t i=0; i<pars.getSize();i++)
	{
		p = (RooRealVar *)pars.at(i);
		outfile<<p->getTitle()<<" "<<p->getVal()<<" "<<p->getError()<<" "<<endl;
	}
	outfile.close();
}

void SysFit::StoreFitResults(string fname, RooFitResult *fitResult)
{
	ofstream outfile;
	outfile.open(fname,std::ios::app);
	time_t t = time(0);   // get time now
    tm* now = localtime(&t);
	//outfile<<" Results of fit performed on: "<<now->tm_mday<<"-"<<(now->tm_mon + 1)<<"-"<<(now->tm_year + 1900)<<"\n"; 
	auto pars = fitResult->floatParsFinal();
	RooRealVar *p;
	for (Int_t i=0; i<pars.getSize();i++)
	{
		p = (RooRealVar *)pars.at(i);
		outfile<<p->getTitle()<<" "<<p->getVal()<<" "<<p->getError()<<" "<<now->tm_mday<<"  "<<(now->tm_mon + 1)<<"  "<<(now->tm_year + 1900)<<" "<<(now->tm_hour)<<":"<<(now->tm_min)<<endl;
	}
	outfile.close();
}

void SysFit::CheckDiscrepancyWrtLastRLcValue(string fname, string channelName)
{
	ifstream fin;
	fin.open(fname);
	string search="RLctau_Isolated";
	string search2;
	//if(GetFitType()=="Simultaneous")
	//	search2="RLcStar";
	//else
	search2="RLcstartau_"+channelName;
	cout<<GetFitType()<<" "<<search2<<endl;
	if(fin)
		cout<<"File Found"<<endl;
	bool isFound=0;
	string temp;
	double RLcTemp=0., sigmaRLcTemp=0.;
	int day=0, month=0, year=0;
	string time;
	vector<double> RLc, sigmaRLc, RLcSt, sigmaRLcSt;
	vector<int> RLcD, RLcM, RLcY, RLcStD, RLcStM, RLcStY;
	vector<string> RLcTime, RLcStTime;
	while(fin>>temp>>RLcTemp>>sigmaRLcTemp>>day>>month>>year>>time)
	{
		if(temp.compare(search)==0)
		{
			RLc.push_back(RLcTemp);
            sigmaRLc.push_back(sigmaRLcTemp);
			RLcD.push_back(day);
			RLcM.push_back(month);
			RLcY.push_back(year);
			RLcTime.push_back(time);
		}
		if(temp.compare(search2)==0)
		{
			RLcSt.push_back(RLcTemp);
            sigmaRLcSt.push_back(sigmaRLcTemp);
			RLcStD.push_back(day);
			RLcStM.push_back(month);
			RLcStY.push_back(year);
			RLcStTime.push_back(time);
		}
	}
	cout<<"-------------- RLc ----------------"<<endl;
	for (int i=0; i<RLc.size();i++)
	{
		cout<<RLc.at(i)<<" "<<sigmaRLc.at(i)<<" "<<RLcD.at(i)<<"-"<<RLcM.at(i)<<"-"<<RLcY.at(i)<<" "<<RLcTime.at(i)<<endl;
	}
	double DeltaRLc=0, DeltaSigmaRLc=0;
	if(RLc.size()>1)
	{
		DeltaRLc = RLc.at(RLc.size()-1) - RLc.at(RLc.size()-2);
		DeltaSigmaRLc = sigmaRLc.at(sigmaRLc.size()-1) - sigmaRLc.at(sigmaRLc.size()-2);
		cout<<"Difference in RLc fitted (blinded) value: "<<DeltaRLc<<" and sigma: " <<DeltaSigmaRLc<<endl;
		cout<<"Compared results of RLc obtained on: "<<RLcD.at(RLc.size()-1)<<"-"<<RLcM.at(RLc.size()-1)<<"-"<<RLcY.at(RLc.size()-1)<<" at "<<RLcTime.at(RLc.size()-1)<<"  and  RLc obtained on: "<<RLcD.at(RLc.size()-2)<<"-"<<RLcM.at(RLc.size()-2)<<"-"<<RLcY.at(RLc.size()-2)<<" at "<<RLcTime.at(RLc.size()-2)<<endl;
	}
	cout<<"-------------- RLcStar ----------------"<<endl;
	for (int i=0; i<RLcSt.size();i++)
	{
		cout<<RLcSt.at(i)<<" "<<sigmaRLcSt.at(i)<<" "<<RLcStD.at(i)<<"-"<<RLcStM.at(i)<<"-"<<RLcStY.at(i)<<" "<<RLcStTime.at(i)<<endl;
	}
	double DeltaRLcSt=0, DeltaSigmaRLcSt=0;
	if(RLcSt.size()>1)
	{
		DeltaRLcSt = RLcSt.at(RLcSt.size()-1) - RLcSt.at(RLcSt.size()-2);
		DeltaSigmaRLcSt = sigmaRLcSt.at(sigmaRLcSt.size()-1) - sigmaRLcSt.at(sigmaRLcSt.size()-2);
		cout<<"Difference in RLcStar fitted (blinded) value: "<<DeltaRLcSt<<" and sigma: " <<DeltaSigmaRLcSt<<endl;
		cout<<"Compared results of RLcStar obtained on: "<<RLcStD.at(RLcSt.size()-1)<<"-"<<RLcStM.at(RLcSt.size()-1)<<"-"<<RLcStY.at(RLcSt.size()-1)<<" at "<<RLcStTime.at(RLcSt.size()-1)<<"  and  RLcStar obtained on: "<<RLcStD.at(RLcSt.size()-2)<<"-"<<RLcStM.at(RLcSt.size()-2)<<"-"<<RLcStY.at(RLcSt.size()-2)<<" at "<<RLcStTime.at(RLcSt.size()-2)<<endl;
	}

}

void SysFit::PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units,string name_suffix, bool legend=kFALSE)
{
	kinemObserv->setRange(title,plotStart,plotEnd); //define the variable you want to plot


	//---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
	RooAbsPdf* model_pdf = model->getPdf(idx->getLabel()); // This gets the pdf in the channel defined by the idx.

	RooArgSet* components_set = model_pdf->getComponents();
	//cout<<"-------- components_set ----------"<<endl;
	// components_set->Print();
	// cout<<"----------------------------------"<<endl;

	TIterator* comp_it = components_set->createIterator();

	RooRealSumPdf* model_sumpdf = new RooRealSumPdf(); // this object will contain the model in a given channel, without the priors constraints

	//--Look for the RooRealSumPdf that contains the model in the given channel
	const char* component_name;
	TObject* comp;
	while ((comp = (TObject*) comp_it->Next()))
	{
		component_name = comp->GetName();
		if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
		{
			cout<<"Component name: "<<component_name<<endl;
			model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
		}
	}

	//Frame for plotting the fit
	RooPlot *frame = new RooPlot();
	frame->SetName(TString(kinemObserv->GetName())+TString("_frame_")+idx->getLabel());
	frame = kinemObserv->frame(RooFit::Title(title));
	//Frame for plotting the pulls
	RooPlot *pframe = new RooPlot();
	pframe->SetName(TString(kinemObserv->GetName())+TString("_pframe_")+idx->getLabel());
	pframe = kinemObserv->frame("");

	//---Take all the components with the name that starts with L_x
	RooArgSet* active_components = model_sumpdf->getComponents();
	std::vector<TString> active_components_names;

	TIterator* active_comp_it  = active_components->createIterator();
	TObject* active_comp;

	while ((active_comp = (TObject*) active_comp_it->Next()))
	{
		TString name(active_comp->GetName());
		if(name.Contains("L_x"))
		{
			active_components_names.push_back(name);
		}
	}

	TString TotComponents = "";

	//Plot the data points
	RooAbsData* channelData = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
	channelData->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));


	//cout<<"******************************************"<<endl;
	for (int i = 0; i < active_components_names.size(); ++i)
	{
		if (TotComponents != "")
		{
			TotComponents+=",";
		}
		TotComponents+="*";
		TotComponents+=active_components_names[i];
		TotComponents+="*";
		Int_t color = GetComponentColor(active_components_names[i]);

		//------- Plot each background component stacking them up
		model_sumpdf->plotOn(frame, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
	}
	//cout<<"******************************************"<<endl;

	//--------Plot the sum of the backgrounds
	model_sumpdf->plotOn(frame, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*data),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));

	//-----Construct the pulls plot
	RooHist *pulls = frame->pullHist();
	pulls->SetFillColor(kGray+1);
	pulls->SetLineColor(kWhite);
	pulls->SetMarkerSize(0.01);

	data->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
	//Build the legend
	if (strcmp(title,"q^{2}")==0)
	{
		Int_t CharmComp=0;
        Int_t CharmStarComp=0;
		TCanvas *clegend = new TCanvas("clegend", "clegend");
		TLegend leg(0.1,0.1,0.9,0.9);
		for (int i=0; i<active_components_names.size(); ++i)
		{
			if(active_components_names[i].Contains("h_2charm-"))
            {
                CharmComp+=1;
                if(CharmComp==1)
                {
					leg.AddEntry((TGraph*)frame->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
                }
            }
			else if(active_components_names[i].Contains("h_starDs-"))
            {
                CharmStarComp+=1;
                if(CharmStarComp==1)
                {
					leg.AddEntry((TGraph*)frame->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
                }
            }
			else
				leg.AddEntry((TGraph*)frame->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
		}
		leg.Draw();
		clegend->Print(TString("plots_")+GetMCcathegory()+TString("/Legend.pdf"));
	}
	
	//-----------Draw the fit result with the pulls plot below
	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->Divide(1,2);

	//----Dimensions of the single pads
	double xmin_1 = 0.005;
	double xmax_1 = 0.995;
	double ymin_1 = 0.35;
	double ymax_1 = 0.995;
	double xmin_2 = .005;
	double xmax_2 = .995;
	double ymin_2 = .05;
	double ymax_2 = .35;

	c1->cd(1)->SetPad(xmin_1, ymin_1, xmax_1, ymax_1 );
	c1->cd(1)->SetBottomMargin(0);
	c1->cd(1)->SetTopMargin(0.15);
	c1->cd(1)->SetLeftMargin(0.14);
	c1->cd(2)->SetPad(xmin_2, ymin_2, xmax_2, ymax_2 );
	c1->cd(2)->SetTopMargin(0);
	c1->cd(2)->SetBottomMargin(0.14/(ymax_2-ymin_2));
	c1->cd(2)->SetLeftMargin(0.14);

	//-draw the bottom pad
	c1->cd(2);
	TPad *padbottom = new TPad("pullspad", "pullspad", 0., 0., 1., 0.3);
	pframe->SetTitle("");

	pframe->GetYaxis()->SetTitle("Pulls");
	pframe->GetYaxis()->SetNdivisions(5);
	pframe->GetXaxis()->SetNdivisions(frame->GetXaxis()->GetNdivisions());
	pframe->GetYaxis()->SetTitleOffset(1.2*(ymax_2-ymin_2));

	pframe->GetYaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe->GetYaxis()->SetLabelSize(0.060/(ymax_2-ymin_2));
	pframe->GetXaxis()->SetLabelSize(0.060/(ymax_2-ymin_2));
	pframe->GetXaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));
	pframe->addPlotable(pulls, "B");

	//pframe->SetMaximum(7.);
	//pframe->SetMinimum(-7.);

	pframe->Draw();

	//--draw the top pad
	c1->cd(1);
	frame->GetYaxis()->SetTitleSize(0.060/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad
	frame->GetYaxis()->SetLabelSize(0.060/(ymax_1-ymin_1));
	frame->GetYaxis()->SetTitleOffset(1.2*(ymax_1-ymin_1));
	frame->SetMinimum(0.01);
	frame->Draw();

	TLatex* lhcb = new TLatex();
	lhcb->SetTextSize(0.06);
	lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");

	c1->Print(TString("plots_")+GetMCcathegory()+TString("/Fit_")+GetFitVarName(title)+TString("_")+TString(name_suffix)+TString("_")+GetFitType()+TString(".pdf"));

}

void SysFit::PlotFitVariables(RooRealVar* fitvar1,const char* title1, RooRealVar* fitvar2, const char* title2, RooRealVar* fitvar3, const char* title3,std::map<std::string,RooDataHist*> data,RooSimultaneous* model, RooCategory* idx,string name_suffix)
{
	fitvar1->SetTitle(title1);
	fitvar2->SetTitle(title2);
	fitvar3->SetTitle(title3);

	//---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
	RooAbsPdf* model_pdf = model->getPdf(idx->getLabel()); // This gets the pdf in the channel defined by the idx.

	RooArgSet* components_set = model_pdf->getComponents();
	cout<<"-------- components_set ----------"<<endl;
	components_set->Print();
	cout<<"----------------------------------"<<endl;

	TIterator* comp_it = components_set->createIterator();

	RooRealSumPdf* model_sumpdf = new RooRealSumPdf(); // this object will contain the model in a given channel, without the priors constraints

	//--Look for the RooRealSumPdf that contains the model in the given channel
	const char* component_name;
	TObject* comp;
	while ((comp = (TObject*) comp_it->Next()))
	{
		component_name = comp->GetName();
		if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
		{
			cout<<"Component name: "<<component_name<<endl;
			model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
		}
	}

	//Frame for plotting the fit variable 1
	RooPlot *frame1 = new RooPlot();
	frame1->SetName(TString(fitvar1->GetName())+TString("_frame_")+idx->getLabel());
	frame1 = fitvar1->frame(RooFit::Title(title1));
	
	//Frame for plotting the fit variable 2
	RooPlot *frame2 = new RooPlot();
	frame2->SetName(TString(fitvar2->GetName())+TString("_frame_")+idx->getLabel());
	frame2 = fitvar2->frame(RooFit::Title(title2));
	
	//Frame for plotting the fit variable 3
	RooPlot *frame3 = new RooPlot();
	frame3->SetName(TString(fitvar3->GetName())+TString("_frame_")+idx->getLabel());
	frame3 = fitvar3->frame(RooFit::Title(title3));
	
	//Frame for plotting the pulls (fitvar1)
	RooPlot *pframe1 = new RooPlot();
	pframe1->SetName(TString(fitvar1->GetName())+TString("_pframe_")+idx->getLabel());
	pframe1 = fitvar1->frame("");
	
	//Frame for plotting the pulls (fitvar2)
	RooPlot *pframe2 = new RooPlot();
	pframe2->SetName(TString(fitvar2->GetName())+TString("_pframe_")+idx->getLabel());
	pframe2 = fitvar2->frame("");

	//Frame for plotting the pulls (fitvar3)
	RooPlot *pframe3 = new RooPlot();
	pframe3->SetName(TString(fitvar3->GetName())+TString("_pframe_")+idx->getLabel());
	pframe3 = fitvar3->frame("");
	
	//---Take all the components with the name that starts with L_x
	RooArgSet* active_components = model_sumpdf->getComponents();
	std::vector<TString> active_components_names;

	TIterator* active_comp_it  = active_components->createIterator();
	TObject* active_comp;

	while ((active_comp = (TObject*) active_comp_it->Next()))
	{
		TString name(active_comp->GetName());
		if(name.Contains("L_x"))
		{
			active_components_names.push_back(name);
		}
	}

	TString TotComponents = "";

	//Plot the data points

	RooAbsData* channelData = (RooAbsData*)data[idx->getLabel()];
	channelData->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));
	channelData->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));
	channelData->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));


	//cout<<"******************************************"<<endl;
	for (int i = 0; i < active_components_names.size(); ++i)
	{
		if (TotComponents != "")
		{
			TotComponents+=",";
		}
		TotComponents+="*";
		TotComponents+=active_components_names[i];
		TotComponents+="*";
		Int_t color = GetComponentColor(active_components_names[i]);

		//------- Plot each background component stacking them up
		model_sumpdf->plotOn(frame1, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
		//------- Plot each background component stacking them up
		model_sumpdf->plotOn(frame2, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
		//------- Plot each background component stacking them up
		model_sumpdf->plotOn(frame3, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
	}
	//cout<<"******************************************"<<endl;

	//--------Plot the sum of the backgrounds
	model_sumpdf->plotOn(frame1, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));
	model_sumpdf->plotOn(frame2, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));
	model_sumpdf->plotOn(frame3, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));

	//-----Construct the pulls plot
	RooHist *pulls1 = frame1->pullHist();
	pulls1->SetFillColor(kGray+1);
	pulls1->SetLineColor(kWhite);
	pulls1->SetMarkerSize(0.01);
	
	RooHist *pulls2 = frame2->pullHist();
	pulls2->SetFillColor(kGray+1);
	pulls2->SetLineColor(kWhite);
	pulls2->SetMarkerSize(0.01);
	
	RooHist *pulls3 = frame3->pullHist();
	pulls3->SetFillColor(kGray+1);
	pulls3->SetLineColor(kWhite);
	pulls3->SetMarkerSize(0.01);

	 data[idx->getLabel()]->plotOn(frame1, RooFit::DrawOption("ZP"));
	 data[idx->getLabel()]->plotOn(frame2, RooFit::DrawOption("ZP"));
	 data[idx->getLabel()]->plotOn(frame3, RooFit::DrawOption("ZP"));


	//data->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
	//data->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
	//data->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));

	
	//-----------Draw the fit result with the pulls plot below
	TCanvas *c1 = new TCanvas("c1", "c1",1500,600);
	c1->Divide(4,2);

	//----Dimensions of the single pads
	double xmin_1[4] = {0.005,0.265,0.525,0.785};
	double xmax_1[4] = {0.26,0.52,0.78,0.995};
	double ymin_1 = 0.35;
	double ymax_1 = 0.995;
	double xmin_2[4] = {0.005,0.265,0.525,0.785};
	double xmax_2[4] = {0.26,0.52,0.78,0.995};
	double ymin_2 = .05;
	double ymax_2 = .35;

	
	for(Int_t i=0;i<4;i++)
	{
		c1->cd(i+1)->SetPad(xmin_1[i], ymin_1, xmax_1[i], ymax_1 );
		c1->cd(i+1)->SetBottomMargin(0);
		c1->cd(i+1)->SetTopMargin(0.15);
		c1->cd(i+1)->SetLeftMargin(0.14);
	}
	for(Int_t i=0;i<4;i++)
	{
		c1->cd(i+5)->SetPad(xmin_2[i], ymin_2, xmax_2[i], ymax_2 );
		c1->cd(i+5)->SetTopMargin(0);
		c1->cd(i+5)->SetBottomMargin(0.14/(ymax_2-ymin_2));
		c1->cd(i+5)->SetLeftMargin(0.14);
	}
	//--draw the top pads
	c1->cd(1);
	frame1 = AdjustVarPlot(frame1, ymin_1, ymax_1, ymin_2, ymax_2);
	frame1->Draw();
	c1->cd(2);
	frame2 = AdjustVarPlot(frame2, ymin_1, ymax_1, ymin_2, ymax_2);
	frame2->Draw();

	c1->cd(3);
	frame3 = AdjustVarPlot(frame3, ymin_1, ymax_1, ymin_2, ymax_2);
	frame3->Draw();
	
	TLatex* lhcb = new TLatex();
	lhcb->SetTextSize(0.04);
	lhcb->DrawLatexNDC(0.6, 0.90, "LHCb Preliminary");

	//-draw the bottom pads
	c1->cd(5);
	TPad *padbottom1 = new TPad("pullspad1", "pullspad1", 0., 0., 1., 0.3);
	pframe1 = AdjustPullsPlot(pframe1,frame1, ymin_1, ymax_1, ymin_2, ymax_2);
	pframe1->addPlotable(pulls1, "B");
	pframe1->Draw();
	
	c1->cd(6);
	TPad *padbottom2 = new TPad("pullspad2", "pullspad2", 0., 0., 1., 0.3);
	pframe2 = AdjustPullsPlot(pframe2, frame2, ymin_1, ymax_1, ymin_2, ymax_2);
	pframe2->addPlotable(pulls2, "B");
	pframe2->Draw();
	c1->cd(7);
	TPad *padbottom3 = new TPad("pullspad3", "pullspad3", 0., 0., 1., 0.3);
	pframe3 = AdjustPullsPlot(pframe3, frame3, ymin_1, ymax_1, ymin_2, ymax_2);
	pframe3->addPlotable(pulls3, "B");
	pframe3->Draw();

	c1->cd(4);
	TPad *padlegend = new TPad("padlegend", "padlegend", 0., 0., 1., 1.);
	gStyle->SetLegendFont(42);
	//Build the legend
	Int_t CharmComp=0;
	Int_t CharmStarComp=0;
	TLegend *leg = new TLegend(0.1,0.05,0.9,0.7);
	for (int i=0; i<active_components_names.size(); ++i)
	{
		/*
		if(active_components_names[i].Contains("h_2charm-"))
		{
			CharmComp+=1;
			if(CharmComp==1)
			{
				leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
			}
		}
		else if(active_components_names[i].Contains("h_starDs-"))
		{
			CharmStarComp+=1;
			if(CharmStarComp==1)
			{
				leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
			}
		}
		else
		*/
			leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
	}
	leg->Draw();


	c1->Print(TString("plots_")+GetMCcathegory()+TString("/FitVariables_")+TString(name_suffix)+TString("_")+GetFitType()+TString(".pdf"));

}

RooPlot* SysFit::AdjustVarPlot(RooPlot *frame, double ymin_1, double ymax_1, double ymin_2, double ymax_2)
{
	frame->GetYaxis()->SetTitleSize(0.040/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad
    frame->GetYaxis()->SetLabelSize(0.025/(ymax_1-ymin_1));
    frame->GetYaxis()->SetTitleOffset(1.6*(ymax_1-ymin_1));
    frame->SetMinimum(0.01);
	return frame;
}

RooPlot* SysFit::AdjustPullsPlot(RooPlot* pframe, RooPlot *frame, double ymin_1, double ymax_1, double ymin_2, double ymax_2)
{
	pframe->SetTitle("");
	pframe->GetYaxis()->SetTitle("Pulls");
	pframe->GetYaxis()->SetNdivisions(5);
	pframe->GetXaxis()->SetNdivisions(frame->GetXaxis()->GetNdivisions());
	pframe->GetYaxis()->SetTitleOffset(1.6*(ymax_2-ymin_2));
	pframe->GetYaxis()->SetTitleSize(0.040/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe->GetYaxis()->SetLabelSize(0.025/(ymax_2-ymin_2));
	pframe->GetXaxis()->SetLabelSize(0.025/(ymax_2-ymin_2));
	pframe->GetXaxis()->SetTitleSize(0.040/(ymax_2-ymin_2));
	return pframe;
}

void SysFit::PlotInBins(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::ModelConfig *mc,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, string name_suffix, bool legend=kFALSE)
{
	//---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
	RooAbsPdf* model_pdf = model->getPdf(idx->getLabel()); // This gets the pdf in the channel defined by the idx.
	RooArgSet* components_set = model_pdf->getComponents();
	TIterator* comp_it = components_set->createIterator();
	//--Look for the RooRealSumPdf that contains the model in the given channel
	RooRealSumPdf* model_sumpdf = new RooRealSumPdf(); // this object will contain the model in a given channel, without the priors constraints
	const char* component_name;
	TObject* comp;
	while ((comp = (TObject*) comp_it->Next()))
	{
		component_name = comp->GetName();
		if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
		{
			//cout<<component_name<<endl;
			model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
		}
	}

	//Frames for plotting the fit
	RooPlot *frame1 = new RooPlot();
	frame1->SetName(TString(kinemObserv->GetName())+TString("_frame1_")+idx->getLabel());
	frame1 = kinemObserv->frame(RooFit::Title(title));
	RooPlot *frame2 = new RooPlot();
	frame2->SetName(TString(kinemObserv->GetName())+TString("_frame2_")+idx->getLabel());
	frame2 = kinemObserv->frame(RooFit::Title(title));
	RooPlot *frame3 = new RooPlot();
	frame3->SetName(TString(kinemObserv->GetName())+TString("_frame3_")+idx->getLabel());
	frame3 = kinemObserv->frame(RooFit::Title(title));
	RooPlot *frame4 = new RooPlot();
	frame4->SetName(TString(kinemObserv->GetName())+TString("_frame4_")+idx->getLabel());
	frame4 = kinemObserv->frame(RooFit::Title(title));

	//Frame for plotting the pulls
	RooPlot *pframe1 = new RooPlot();
	pframe1->SetName(TString(kinemObserv->GetName())+TString("_pframe1_")+idx->getLabel());
	pframe1 = kinemObserv->frame("");
	RooPlot *pframe2 = new RooPlot();
	pframe2->SetName(TString(kinemObserv->GetName())+TString("_pframe2_")+idx->getLabel());
	pframe2 = kinemObserv->frame("");
	RooPlot *pframe3 = new RooPlot();
	pframe3->SetName(TString(kinemObserv->GetName())+TString("_pframe3_")+idx->getLabel());
	pframe3 = kinemObserv->frame("");
	RooPlot *pframe4 = new RooPlot();
	pframe4->SetName(TString(kinemObserv->GetName())+TString("_pframe4_")+idx->getLabel());
	pframe4 = kinemObserv->frame("");

	//---Take all the components with the name that starts with L_x
	RooArgSet* active_components = model_sumpdf->getComponents();
	std::vector<TString> active_components_names;

	TIterator* active_comp_it  = active_components->createIterator();
	TObject* active_comp;

	while ((active_comp = (TObject*) active_comp_it->Next()))
	{
		TString name(active_comp->GetName());
		//cout<<name<<endl;
		if(name.Contains("L_x"))
		{
			active_components_names.push_back(name);
		}
	}

	// Get the variable the projection is defined on
	RooRealVar* q2 = (RooRealVar*)mc->GetObservables()->find((std::string("obs_x_RLc_kinematic_")+std::string(name_suffix)).c_str());
	q2->setRange("q2_range1", -2.,    2.);
	q2->setRange("q2_range2", 2., 6.);
	q2->setRange("q2_range3", 6.,  10.);
	q2->setRange("q2_range4", 10., 14.);

	//Plot data points
	RooAbsData* channelData1 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
	RooAbsData* channelData2 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
	RooAbsData* channelData3 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
	RooAbsData* channelData4 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));

	channelData1->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range1"));
	channelData2->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range2"));
	channelData3->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range3"));
	channelData4->plotOn(frame4, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range4"));

	//Plot all components
	TString TotComponents = "";
	for (int i = 0; i < active_components_names.size(); ++i)
	{
		TString component = "";
		component+="*";
		component+=active_components_names[i];
		component+="*";
		TotComponents+=component;
		if(TotComponents!="")
			TotComponents+=",";

		Int_t color = GetComponentColor(active_components_names[i]);
		//------- Plot each background component stacking them up
		model_sumpdf->plotOn(frame1, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData1),RooFit::ProjectionRange("q2_range1"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot1").Data()));
		model_sumpdf->plotOn(frame2, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData2),RooFit::ProjectionRange("q2_range2"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot2").Data()));
		model_sumpdf->plotOn(frame3, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData3),RooFit::ProjectionRange("q2_range3"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot3").Data()));
		model_sumpdf->plotOn(frame4, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData4),RooFit::ProjectionRange("q2_range4"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot4").Data()));


	}

	//--------Plot the sum of the backgrounds
	model_sumpdf->plotOn(frame1, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData1),RooFit::ProjectionRange("q2_range1"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents));
	model_sumpdf->plotOn(frame2, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData2),RooFit::ProjectionRange("q2_range2"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents));
	model_sumpdf->plotOn(frame3, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData3),RooFit::ProjectionRange("q2_range3"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents));
	model_sumpdf->plotOn(frame4, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData4),RooFit::ProjectionRange("q2_range4"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents));


	channelData1->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range1"));
	channelData2->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range2"));
	channelData3->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range3"));
	channelData4->plotOn(frame4, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson),RooFit::CutRange("q2_range4"));

	//-----Construct the pulls plot


	RooHist *pulls1 = frame1->pullHist();
	RooHist *pulls2 = frame2->pullHist();
	RooHist *pulls3 = frame3->pullHist();
	RooHist *pulls4 = frame4->pullHist();

	pulls1->SetFillColor(kGray+1);
	pulls2->SetFillColor(kGray+1);
	pulls3->SetFillColor(kGray+1);
	pulls4->SetFillColor(kGray+1);

	pulls1->SetLineColor(kWhite);
	pulls2->SetLineColor(kWhite);
	pulls3->SetLineColor(kWhite);
	pulls4->SetLineColor(kWhite);

	pulls1->SetMarkerSize(0.01);
	pulls2->SetMarkerSize(0.01);
	pulls3->SetMarkerSize(0.01);
	pulls4->SetMarkerSize(0.01);

	pframe1->addPlotable(pulls1, "B");
	pframe2->addPlotable(pulls2, "B");
	pframe3->addPlotable(pulls3, "B");
	pframe4->addPlotable(pulls4, "B");

	//-----------Draw the fit result with the pulls plot below

	TCanvas *c1 = new TCanvas("c1", "c1");
	TCanvas *c2 = new TCanvas("c2", "c2");
	TCanvas *c3 = new TCanvas("c3", "c3");
	TCanvas *c4 = new TCanvas("c4", "c4");

	c1->Divide(1,2); 
	c2->Divide(1,2); 
	c3->Divide(1,2); 
	c4->Divide(1,2); 

	//----Dimensions of the single pads
	double xmin_1 = 0.005;
	double xmax_1 = 0.995;
	double ymin_1 = 0.35;
	double ymax_1 = 0.995;
	double xmin_2 = .005;
	double xmax_2 = .995;
	double ymin_2 = .05;  
	double ymax_2 = .35;

	c1->cd(1)->SetPad(xmin_1, ymin_1, xmax_1, ymax_1 );
	c2->cd(1)->SetPad(xmin_1, ymin_1, xmax_1, ymax_1 );
	c3->cd(1)->SetPad(xmin_1, ymin_1, xmax_1, ymax_1 );
	c4->cd(1)->SetPad(xmin_1, ymin_1, xmax_1, ymax_1 );

	c1->cd(1)->SetBottomMargin(0);
	c2->cd(1)->SetBottomMargin(0);
	c3->cd(1)->SetBottomMargin(0);
	c4->cd(1)->SetBottomMargin(0);

	c1->cd(1)->SetTopMargin(0.15);
	c2->cd(1)->SetTopMargin(0.15);
	c3->cd(1)->SetTopMargin(0.15);
	c4->cd(1)->SetTopMargin(0.15);

	c1->cd(1)->SetLeftMargin(0.14);
	c2->cd(1)->SetLeftMargin(0.14);
	c3->cd(1)->SetLeftMargin(0.14);
	c4->cd(1)->SetLeftMargin(0.14);

	c1->cd(2)->SetPad(xmin_2, ymin_2, xmax_2, ymax_2 );
	c2->cd(2)->SetPad(xmin_2, ymin_2, xmax_2, ymax_2 );
	c3->cd(2)->SetPad(xmin_2, ymin_2, xmax_2, ymax_2 );
	c4->cd(2)->SetPad(xmin_2, ymin_2, xmax_2, ymax_2 );

	c1->cd(2)->SetTopMargin(0);
	c2->cd(2)->SetTopMargin(0);
	c3->cd(2)->SetTopMargin(0);
	c4->cd(2)->SetTopMargin(0);

	c1->cd(2)->SetBottomMargin(0.14/(ymax_2-ymin_2));
	c2->cd(2)->SetBottomMargin(0.14/(ymax_2-ymin_2));
	c3->cd(2)->SetBottomMargin(0.14/(ymax_2-ymin_2));
	c4->cd(2)->SetBottomMargin(0.14/(ymax_2-ymin_2));

	c1->cd(2)->SetLeftMargin(0.14);
	c2->cd(2)->SetLeftMargin(0.14);
	c3->cd(2)->SetLeftMargin(0.14);
	c4->cd(2)->SetLeftMargin(0.14);

	//-draw the bottom pad 
	TPad *padbottom = new TPad("pullspad", "pullspad", 0., 0., 1., 0.3);

	pframe1->SetTitle("");
	pframe2->SetTitle("");
	pframe3->SetTitle("");
	pframe4->SetTitle("");

	pframe1->GetYaxis()->SetTitle("Pulls");
	pframe2->GetYaxis()->SetTitle("Pulls");
	pframe3->GetYaxis()->SetTitle("Pulls");
	pframe4->GetYaxis()->SetTitle("Pulls");

	pframe1->GetYaxis()->SetNdivisions(5);
	pframe2->GetYaxis()->SetNdivisions(5);
	pframe3->GetYaxis()->SetNdivisions(5);
	pframe4->GetYaxis()->SetNdivisions(5);

	pframe1->GetXaxis()->SetNdivisions(frame1->GetXaxis()->GetNdivisions());
	pframe2->GetXaxis()->SetNdivisions(frame2->GetXaxis()->GetNdivisions());
	pframe3->GetXaxis()->SetNdivisions(frame3->GetXaxis()->GetNdivisions());
	pframe4->GetXaxis()->SetNdivisions(frame4->GetXaxis()->GetNdivisions());

	pframe1->GetYaxis()->SetTitleOffset(1.2*(ymax_2-ymin_2));
	pframe2->GetYaxis()->SetTitleOffset(1.2*(ymax_2-ymin_2));
	pframe3->GetYaxis()->SetTitleOffset(1.2*(ymax_2-ymin_2));
	pframe4->GetYaxis()->SetTitleOffset(1.2*(ymax_2-ymin_2));

	pframe1->GetYaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe1->GetXaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe2->GetYaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe2->GetXaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe3->GetYaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe3->GetXaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe4->GetYaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad
	pframe4->GetXaxis()->SetTitleSize(0.060/(ymax_2-ymin_2));//rescale the title and the label sizes after dividing the pad

	pframe1->GetYaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe1->GetXaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe2->GetYaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe2->GetXaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe3->GetYaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe3->GetXaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe4->GetYaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 
	pframe4->GetXaxis()->SetLabelSize(0.060/(ymax_2-ymin_2)); 

	pframe1->SetMaximum(7.);  
	pframe2->SetMaximum(7.);  
	pframe3->SetMaximum(7.);  
	pframe4->SetMaximum(7.);  

	pframe1->SetMinimum(-7.);
	pframe2->SetMinimum(-7.);
	pframe3->SetMinimum(-7.);
	pframe4->SetMinimum(-7.);

	//Plot the pulls on each bottom canvas
	c1->cd(2);
	pframe1->Draw();

	c2->cd(2);
	pframe2->Draw();

	c3->cd(2);
	pframe3->Draw();

	c4->cd(2);
	pframe4->Draw();

	//--draw the top pad
	frame1->GetYaxis()->SetTitleSize(0.060/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad
	frame2->GetYaxis()->SetTitleSize(0.060/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad
	frame3->GetYaxis()->SetTitleSize(0.060/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad
	frame4->GetYaxis()->SetTitleSize(0.060/(ymax_1-ymin_1)); //rescale the title and the label sizes after dividing the pad

	frame1->GetYaxis()->SetLabelSize(0.060/(ymax_1-ymin_1));
	frame2->GetYaxis()->SetLabelSize(0.060/(ymax_1-ymin_1));
	frame3->GetYaxis()->SetLabelSize(0.060/(ymax_1-ymin_1));
	frame4->GetYaxis()->SetLabelSize(0.060/(ymax_1-ymin_1));

	frame1->GetYaxis()->SetTitleOffset(1.2*(ymax_1-ymin_1));
	frame2->GetYaxis()->SetTitleOffset(1.2*(ymax_1-ymin_1));
	frame3->GetYaxis()->SetTitleOffset(1.2*(ymax_1-ymin_1));
	frame4->GetYaxis()->SetTitleOffset(1.2*(ymax_1-ymin_1));

	frame1->SetMinimum(0.01);
	frame2->SetMinimum(0.01);
	frame3->SetMinimum(0.01);
	frame4->SetMinimum(0.01);

	c1->cd(1);
	frame1->Draw();
	c2->cd(1);
	frame2->Draw();
	c3->cd(1);
	frame3->Draw();
	c4->cd(1);
	frame4->Draw();

	TLatex* lhcb = new TLatex();
	lhcb->SetTextSize(0.06);

	TLatex* q2bin1 = new TLatex();
	TLatex* q2bin2 = new TLatex();
	TLatex* q2bin3 = new TLatex();
	TLatex* q2bin4 = new TLatex();

	c1->cd(1);
	lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
	q2bin1->DrawLatexNDC(0.4, 0.90, "-2. < q^{2} < 2. GeV^{2}/c^{4}");
	c2->cd(1);
	lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
	q2bin2->DrawLatexNDC(0.4, 0.90, "2. < q^{2} < 6. GeV^{2}/c^{4}");
	c3->cd(1);
	lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
	q2bin3->DrawLatexNDC(0.4, 0.90, "6. < q^{2} < 10. GeV^{2}/c^{4}");
	c4->cd(1);
	lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
	q2bin4->DrawLatexNDC(0.4, 0.90, "10. < q^{2} < 14. GeV^{2}/c^{4}");

	//Save canvases
	c1->Print(TString("plots_q2bin/Fit_")+GetFitVarName(title)+TString("_q2bin1_")+TString(name_suffix)+TString(".pdf"));
	c2->Print(TString("plots_q2bin/Fit_")+GetFitVarName(title)+TString("_q2bin2_")+TString(name_suffix)+TString(".pdf"));
	c3->Print(TString("plots_q2bin/Fit_")+GetFitVarName(title)+TString("_q2bin3_")+TString(name_suffix)+TString(".pdf"));
	c4->Print(TString("plots_q2bin/Fit_")+GetFitVarName(title)+TString("_q2bin4_")+TString(name_suffix)+TString(".pdf"));
}

void SysFit::PlotFitVariablesInFit(RooRealVar* fitvar1,const char* title1, RooRealVar* fitvar2, const char* title2, RooRealVar* fitvar3, const char* title3,RooAbsData* data,RooSimultaneous* model, RooCategory* idx,string name_suffix)
{
    fitvar1->SetTitle(title1);
    fitvar2->SetTitle(title2);
    fitvar3->SetTitle(title3);

    //---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
    RooAbsPdf* model_pdf = model->getPdf(idx->getLabel()); // This gets the pdf in the channel defined by the idx.

    RooArgSet* components_set = model_pdf->getComponents();
    //cout<<"-------- components_set ----------"<<endl;
    // components_set->Print();
    // cout<<"----------------------------------"<<endl;

    TIterator* comp_it = components_set->createIterator();

    RooRealSumPdf* model_sumpdf = new RooRealSumPdf(); // this object will contain the model in a given channel, without the priors constraints

    //--Look for the RooRealSumPdf that contains the model in the given channel
    const char* component_name;
    TObject* comp;
    while ((comp = (TObject*) comp_it->Next()))
    {
        component_name = comp->GetName();
        if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
        {
            cout<<"Component name: "<<component_name<<endl;
            model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
        }
    }

    //Frame for plotting the fit variable 1
    RooPlot *frame1 = new RooPlot();
    frame1->SetName(TString(fitvar1->GetName())+TString("_frame_")+idx->getLabel());
    frame1 = fitvar1->frame(RooFit::Title(title1));

    //Frame for plotting the fit variable 2
    RooPlot *frame2 = new RooPlot();
    frame2->SetName(TString(fitvar2->GetName())+TString("_frame_")+idx->getLabel());
    frame2 = fitvar2->frame(RooFit::Title(title2));

    //Frame for plotting the fit variable 3
	RooPlot *frame3 = new RooPlot();
    frame3->SetName(TString(fitvar3->GetName())+TString("_frame_")+idx->getLabel());
    frame3 = fitvar3->frame(RooFit::Title(title3));

    //Frame for plotting the pulls (fitvar1)
    RooPlot *pframe1 = new RooPlot();
    pframe1->SetName(TString(fitvar1->GetName())+TString("_pframe_")+idx->getLabel());
    pframe1 = fitvar1->frame("");

    //Frame for plotting the pulls (fitvar2)
    RooPlot *pframe2 = new RooPlot();
    pframe2->SetName(TString(fitvar2->GetName())+TString("_pframe_")+idx->getLabel());
    pframe2 = fitvar2->frame("");

    //Frame for plotting the pulls (fitvar3)
    RooPlot *pframe3 = new RooPlot();
    pframe3->SetName(TString(fitvar3->GetName())+TString("_pframe_")+idx->getLabel());
    pframe3 = fitvar3->frame("");

    //---Take all the components with the name that starts with L_x
    RooArgSet* active_components = model_sumpdf->getComponents();
    std::vector<TString> active_components_names;

    TIterator* active_comp_it  = active_components->createIterator();
    TObject* active_comp;

    while ((active_comp = (TObject*) active_comp_it->Next()))
    {
        TString name(active_comp->GetName());
        if(name.Contains("L_x"))
        {
            active_components_names.push_back(name);
        }
    }

    TString TotComponents = "";
	//Plot the data points
    RooAbsData* channelData = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
    channelData->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));
    channelData->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));
    channelData->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));


    //cout<<"******************************************"<<endl;
    for (int i = 0; i < active_components_names.size(); ++i)
    {
        if (TotComponents != "")
        {
            TotComponents+=",";
        }
        TotComponents+="*";
        TotComponents+=active_components_names[i];
        TotComponents+="*";
        Int_t color = GetComponentColor(active_components_names[i]);

        //------- Plot each background component stacking them up
        model_sumpdf->plotOn(frame1, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
        //------- Plot each background component stacking them up
        model_sumpdf->plotOn(frame2, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
        //------- Plot each background component stacking them up
        model_sumpdf->plotOn(frame3, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot")));
    }
    //cout<<"******************************************"<<endl;

    //--------Plot the sum of the backgrounds
    model_sumpdf->plotOn(frame1, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*data),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));
    model_sumpdf->plotOn(frame2, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*data),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));
    model_sumpdf->plotOn(frame3, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*data),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kBlack));

	//-----Construct the pulls plot
    RooHist *pulls1 = frame1->pullHist();
    pulls1->SetFillColor(kGray+1);
    pulls1->SetLineColor(kWhite);
    pulls1->SetMarkerSize(0.01);

    RooHist *pulls2 = frame2->pullHist();
    pulls2->SetFillColor(kGray+1);
    pulls2->SetLineColor(kWhite);
    pulls2->SetMarkerSize(0.01);

    RooHist *pulls3 = frame3->pullHist();
    pulls3->SetFillColor(kGray+1);
    pulls3->SetLineColor(kWhite);
    pulls3->SetMarkerSize(0.01);

    data->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
    data->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
    data->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));


    //-----------Draw the fit result with the pulls plot below
    TCanvas *c1 = new TCanvas("c1", "c1",1500,600);
    c1->Divide(4,2);

    //----Dimensions of the single pads
    double xmin_1[4] = {0.005,0.265,0.525,0.785};
    double xmax_1[4] = {0.26,0.52,0.78,0.995};
    double ymin_1 = 0.35;
    double ymax_1 = 0.995;
    double xmin_2[4] = {0.005,0.265,0.525,0.785};
	double xmax_2[4] = {0.26,0.52,0.78,0.995};
	double ymin_2 = .05;
    double ymax_2 = .35;


    for(Int_t i=0;i<4;i++)
    {
        c1->cd(i+1)->SetPad(xmin_1[i], ymin_1, xmax_1[i], ymax_1 );
        c1->cd(i+1)->SetBottomMargin(0);
        c1->cd(i+1)->SetTopMargin(0.15);
        c1->cd(i+1)->SetLeftMargin(0.14);
    }
    for(Int_t i=0;i<4;i++)
    {
        c1->cd(i+5)->SetPad(xmin_2[i], ymin_2, xmax_2[i], ymax_2 );
        c1->cd(i+5)->SetTopMargin(0);
        c1->cd(i+5)->SetBottomMargin(0.14/(ymax_2-ymin_2));
        c1->cd(i+5)->SetLeftMargin(0.14);
    }
    //--draw the top pads
    c1->cd(1);
    frame1 = AdjustVarPlot(frame1, ymin_1, ymax_1, ymin_2, ymax_2);
    frame1->Draw();
    c1->cd(2);
    frame2 = AdjustVarPlot(frame2, ymin_1, ymax_1, ymin_2, ymax_2);
    frame2->Draw();

    c1->cd(3);
    frame3 = AdjustVarPlot(frame3, ymin_1, ymax_1, ymin_2, ymax_2);
    frame3->Draw();

	 TLatex* lhcb = new TLatex();
    lhcb->SetTextSize(0.04);
    lhcb->DrawLatexNDC(0.6, 0.90, "LHCb Preliminary");

    //-draw the bottom pads
    c1->cd(5);
    TPad *padbottom1 = new TPad("pullspad1", "pullspad1", 0., 0., 1., 0.3);
    pframe1 = AdjustPullsPlot(pframe1,frame1, ymin_1, ymax_1, ymin_2, ymax_2);
    pframe1->addPlotable(pulls1, "B");
    pframe1->Draw();

    c1->cd(6);
    TPad *padbottom2 = new TPad("pullspad2", "pullspad2", 0., 0., 1., 0.3);
    pframe2 = AdjustPullsPlot(pframe2, frame2, ymin_1, ymax_1, ymin_2, ymax_2);
    pframe2->addPlotable(pulls2, "B");
    pframe2->Draw();

    c1->cd(7);
    TPad *padbottom3 = new TPad("pullspad3", "pullspad3", 0., 0., 1., 0.3);
    pframe3 = AdjustPullsPlot(pframe3, frame3, ymin_1, ymax_1, ymin_2, ymax_2);
    pframe3->addPlotable(pulls3, "B");
    pframe3->Draw();

    c1->cd(4);
    TPad *padlegend = new TPad("padlegend", "padlegend", 0., 0., 1., 1.);
    gStyle->SetLegendFont(42);
    //Build the legend
    Int_t CharmComp=0;
    Int_t CharmStarComp=0;
    TLegend *leg = new TLegend(0.1,0.05,0.9,0.7);
	for (int i=0; i<active_components_names.size(); ++i)
    {
        if(active_components_names[i].Contains("h_2charm-"))
        {
            CharmComp+=1;
            if(CharmComp==1)
            {
                leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
            }
        }
        else if(active_components_names[i].Contains("h_starDs-"))
        {
            CharmStarComp+=1;
            if(CharmStarComp==1)
            {
                leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
            }
        }
        else
            leg->AddEntry((TGraph*)frame1->findObject(active_components_names[i]+TString("_plot").Data()),GetComponentName(active_components_names[i]),"F");
    }
    leg->Draw();


    c1->Print(TString("plots_")+GetMCcathegory()+TString("/FitVariables_direct_")+TString(name_suffix)+TString("_")+GetFitType()+TString(".pdf"));

}



