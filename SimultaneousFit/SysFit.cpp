#include "SysFit.h"
#include <TStyle.h>
#include <TRandom3.h>

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
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncmu_Isolated",{8.E5,3.E3,2.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nc2charm-2body_Isolated",{4.0E4,3.E3,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nc2charm-mbody_Isolated",{4.0E4,3.E3,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncstarmu_Isolated",{2.6E5,3.E3,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcstarDs_Isolated",{2E4,3.E3,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcMISID_Isolated",{3.1E4,10,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("NcCombinatorial_Isolated",{4.0E4,1.E1,1.E6}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Nctau_Isolated",{0.05,1.E-9,1.}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncstartau_Isolated",{0.05,1.E-9,1.}));
	start_parameters["Isolated"].insert(pair<string, vector<Double_t>>("Ncpha_Isolated",{0,-3.,3.}));
	
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncmu_Kenriched",{1.4E3,100,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncstarmu_Kenriched",{1.E3,10,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcstarDs_Kenriched",{1.E3,10,1.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Ncstartau_Kenriched",{0.,1E-9,1}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nctau_Kenriched",{0.,1E-9,1.}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nc2charm-2body_Kenriched",{5.E3,500,2.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("Nc2charm-mbody_Kenriched",{5.E3,500,2.E4}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcCombinatorial_Kenriched",{3.2E3,1.E2,5.E3}));
	start_parameters["Kenriched"].insert(pair<string, vector<Double_t>>("NcMISID_Kenriched",{4.6E3,100,2.E4}));
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

/*
void SysFit::SetStartParameters(map<string,vector<Double_t> > parameters)
{
	for(std::map<std::string,std::vector<Double_t> > ::iterator iP= parameters.begin(); iP != parameters.end();++iP)
	{
		start_parameters[iP->first.substr(2)][0] = iP->second[0];
	}
}*/

void SysFit::SetStartParameters(RooFitResult *fitResult, string ch_name)
{
	RooArgList ParList = fitResult->floatParsFinal();

	for(Int_t i=0;i<ParList.getSize();i++)
	{
		const char* title =  ParList.at(i)->GetTitle();
		if(string(title)!="RLctau_"+ch_name)
		{
		RooRealVar* value = (RooRealVar*)(ParList.find(title));
		RooErrorVar* error = value->errorVar();
		//start_parameters[ch_name].insert(pair<string, vector<Double_t>>(string(title),{value->getVal(), value->getVal()-error->getVal(),value->getVal()+error->getVal() }));
		start_parameters[ch_name][string(title)] = {value->getVal(), value->getVal()-3*error->getVal(),value->getVal()+3*error->getVal() };
		}
	}
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
	if (component.Contains("h_mu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #mu #nu";
	else if(component.Contains("h_starmu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}^{*} #mu #nu_{#mu}";
	else if(component.Contains("h_starDs_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}^{*} X_{c}";
	else if(component.Contains("h_startau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}^{*} #tau #nu_{#tau}";
	else if(component.Contains("h_tau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #tau #nu_{#tau}";
	else if(component.Contains("h_2charm")) name = "#Lambda_{b} #rightarrow #Lambda_{c} X_{c}";
	else if(component.Contains("h_MISID")) name = "MISID";
	else if(component.Contains("h_Combinatorial")) name =  "Combinatorial";

	return name;
}

Int_t SysFit::GetComponentColor(TString component)
{       
	Int_t color=-99;
	if (component.Contains("h_mu_")) color = kBlue;
	else if(component.Contains("h_starmu_")) color = kViolet;
	else if(component.Contains("h_startau_")) color = kMagenta;
	else if(component.Contains("h_starDs_")) color = kGreen+10;
	else if(component.Contains("h_tau_"))color = kRed;
	else if(component.Contains("h_2charm")) color = kGreen;
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
	//Read root file with histograms of samples
	TFile histFile(inputFile.c_str());
	TH1 *h_temp;
	//Get the histo from the file and save in htemp
	histFile.GetObject(("h_"+type).c_str(),h_temp);
	//Verify the histo is non null
	assert(h_temp!=NULL);

	//Get Normalisation factor to 1
	Double_t mcNorm = 1./h_temp->Integral();

	cout<<" Integral_"+type+" = "<<h_temp->Integral()<<endl;
	cout<<" Normalisation factor "+type+" = " << mcNorm <<endl;

	return mcNorm;
}



void SysFit::AddSample(string type, string inputFile, bool shapeUncert, bool GaussConstraints, const bool BBeast, Channel** chan, vector<Double_t> params, vector<Double_t> mu_params={0})
{
	//Define the sample
	Sample sample;
	//set the file from where to take the template
	sample.SetInputFile(inputFile);
	//get the normalisation factor to 1
	Double_t mcNorms = GetHistoNormalisation(inputFile, type);
	//Get the name of the channel
	string chname =(*chan)->GetName();
	string category;
	if (chname.find("Isolated") != std::string::npos) category = "Isolated";
	if (chname.find("Kenriched") != std::string::npos) category = "Kenriched";


	//Code for shape variations for 2charm sample
	if(shapeUncert && type.find("2charm-mbody") !=std::string::npos)
	{
		sample.SetName("h_"+type);
		sample.SetHistoName("h_"+type);
		sample.AddHistoSys(type+"_variation","h_"+type+"_1ml",inputFile,"", "h_"+type+"_1pl",inputFile,"");

		sample.AddHistoSys(type+"_quadratic_variation","h_"+type+"_1mq",inputFile,"", "h_"+type+"_1pq",inputFile,"");
	}
	if(shapeUncert && type.find("2charm") == std::string::npos)
	{
		sample.SetName("h_w_"+type+"_mean");
		sample.SetHistoName("h_w_"+type+"_mean");
		sample.AddHistoSys(type+"_shape_unc","h_w_"+type+"_min",inputFile,"","h_w_"+type+"_max",inputFile,"");
	}
	if(!shapeUncert)
	{
		sample.SetName("h_"+type);
		sample.SetHistoName("h_"+type);
	}
	//set a sample to be "normalised by theory" (its normalisation scales with luminosity)
	sample.SetNormalizeByTheory(kFALSE);

	//add the normalisation factors which will be multiplied
	if(type=="tau" || type=="startau") 
	{
		sample.AddNormFactor("Ncmu_"+category,mu_params[0],mu_params[1],mu_params[2]);
		sample.AddNormFactor("RLc"+type+"_"+category, params[0],params[1],params[2]);
	}
	else
		sample.AddNormFactor("Nc"+type+"_"+category, params[0],params[1],params[2]);
	sample.AddNormFactor("mcNc"+type+"_"+category,mcNorms,1e-20,1);

	//----- Add Gaussian Constraint for MISID/Combinatorial for control fit part
	if (GaussConstraints && chname=="RLc_kinematic_Kenriched")
	{
		Double_t w = GetWeightGaussConstraint(type);
		sample.AddOverallSys("NormUnc_"+type+"_"+category,1.,w);
	}

	if(BBeast) sample.ActivateStatError();
	(*chan)->AddSample(sample);
}

void SysFit::PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooSimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units,string name_suffix, bool legend=kFALSE)
{
        kinemObserv->setRange(title,plotStart,plotEnd); //define the variable you want to plot


        //---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
        RooAbsPdf* model_pdf = model->getPdf(idx->getLabel()); // This gets the pdf in the channel defined by the idx.

        RooArgSet* components_set = model_pdf->getComponents();
        //cout<<"-------- components_set ----------"<<endl;
        //components_set->Print();
        //cout<<"----------------------------------"<<endl;
        
		TIterator* comp_it = components_set->createIterator();

        RooRealSumPdf* model_sumpdf = new RooRealSumPdf(); // this object will contain the model in a given channel, without the priors constraints

        //--Look for the RooRealSumPdf that contains the model in the given channel
        const char* component_name;
        TObject* comp;
	while (comp = (TObject*) comp_it->Next())
        {
                component_name = comp->GetName();
                if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
                {
                        //cout<<component_name<<endl;
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

        while (active_comp = (TObject*) active_comp_it->Next())
        {
                TString name(active_comp->GetName());
                //cout<<name<<endl;
                if(name.Contains("L_x"))
                {
                        active_components_names.push_back(name);
                }
        }

        TString TotComponents = "";

	//Plot the data points
        //data->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
	RooAbsData* channelData = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  channelData->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError(RooAbsData::Poisson));


        //cout<<"******************************************"<<endl;
        for (int i = 0; i < active_components_names.size(); ++i)
        {
                //std::cout << active_components_names[i] << std::endl;
                TString component = "";
                component+="*";
                component+=active_components_names[i];
                component+="*";
                TotComponents+=component;
                if(TotComponents!="")
                        TotComponents+=",";

                Int_t color = GetComponentColor(active_components_names[i]);
                //cout<<color<<endl;
                //cout<<component<<"     "<<TotComponents<<endl;
                //
                //------- Plot each background component
                //model_sumpdf->plotOn(frame, RooFit::Slice(*idx),RooFit::ProjWData(*idx,*data),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(component),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));
                //------- Plot each background component stacking them up
                model_sumpdf->plotOn(frame, RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));
        }
        //cout<<"******************************************"<<endl;

        //--------Plot the sum of the backgrounds
        model_sumpdf->plotOn(frame, RooFit::Slice(*idx,TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents)), RooFit::MoveToBack(),RooFit::Range("title");

data->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
        //-----Construct the pulls plot
        RooHist *pulls = frame->pullHist();
        pulls->SetFillColor(kGray+1);
	//Build the legend
        TCanvas *clegend = new TCanvas("clegend", "clegend");
        TLegend leg(0.1,0.1,0.9,0.9);

        if (strcmp(title,"q^{2}")==0)
        {
                std::vector<TGraph*> component_graphs;
                Int_t CharmComp=0;
                for (int i=0; i<active_components_names.size(); ++i)
                {
                        //cout<<endl;
                        //cout<<active_components_names[i]<<endl;
                        Int_t c = frame->numItems() -3 -i;
                        //cout<<c<<endl;
                        //frame->getObject(c)->Print();
                        //cout<<frame->numItems()<<endl;
                        //      cout<<endl;

                        TGraph* graph = (TGraph*) frame->getObject(frame->numItems()-3-i);
                        if(active_components_names[i].Contains("h_2charm_"))
                        {
                                CharmComp+=1;
                                if(CharmComp==1)
                                        leg.AddEntry(graph, GetComponentName(active_components_names[i]), "f");
                                else
                                        continue;
                        }
                        else
                                 leg.AddEntry(graph, GetComponentName(active_components_names[i]), "f");
                }
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
        if (strcmp(title,"q^{2}")==0)
        {
                clegend->cd();
                leg.Draw();
                clegend->Print("plots/Legend.pdf");
                //clegend->Print("plots/Legend.C");
        }

        c1->Print(TString("plots/Fit_")+GetFitVarName(title)+TString("_")+TString(name_suffix)+TString(".pdf"));

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
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcmu_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNctau_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarmu_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstartau_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarDs_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm-mbody_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm-2body_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcMISID_"+channel).c_str())))->setConstant(kTRUE);
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcCombinatorial_"+channel).c_str())))->setConstant(kTRUE);
	
	return mc;
}

RooStats::ModelConfig* SysFit::FixYields(RooStats::ModelConfig *mc, string sample, string channel)
{
	((RooRealVar*)(mc->GetNuisanceParameters()->find(("Nc"+sample+"_"+channel).c_str())))->setConstant(kTRUE);
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
	Int_t nchannels = channel_names.size();
	vector <Channel*> channels;
	for(Int_t i=0; i<nchannels; i++)
	{
		channels.push_back(new Channel((string("RLc_kinematic_")+channel_names[i]).c_str()));
		channels[i]->SetStatErrorConfig(1e-5, "Poisson");
		map<string,vector<Double_t>> start_param = GetStartParameters(channel_names[i]);
		vector<string> category = GetCategory(start_param);
		vector<string> param_names = GetParametersName(start_param);
		string filename = string("RootFiles/Histos_")+channel_names[i]+string(".root");

		//Add the samples to the channels
		for(Int_t j=0; j<category.size(); j++)
		{
			cout<<category[j]<<endl;

			if(category[j]=="pha")
				continue;
			if(category[j]!="tau")
				AddSample(category[j],filename,IsShapeUncertain(category[j]),IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]]);
			if(category[j]=="tau")//NB:: To modify the 0 adding ShapeUnc
				AddSample(category[j],filename,0,IsGaussConstrained(category[j]),BBeast,&channels[i],start_param[param_names[j]],start_param[string("Ncmu_")+channel_names[i]]);
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
	w->Print();
	cout<<"-----------------------------------------------------------"<<endl;
	return w;
}


RooStats::ModelConfig* SysFit::CreateModel(RooWorkspace* w)
{
	//Get the model manually
	RooStats::ModelConfig *mc = (RooStats::ModelConfig*) w->obj("ModelConfig");
	return mc;
}

RooFitResult* SysFit::Fit(RooStats::ModelConfig* mc, RooStats::HistFactory::Measurement meas, RooWorkspace *w)
{
	RooSimultaneous *model_ = (RooSimultaneous*)mc->GetPdf();

	//Fix the MC normalisation of the histos
	cout << "-------------------------   Fixing parameters   ------------------------------" << endl;
	//mc->GetNuisanceParameters()->Print();
	((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);
	vector <string> channel_names = NameChannels();
	Int_t nchannels = channel_names.size();
	for(Int_t i=0; i<nchannels; i++)
		SetChannelConstants(mc,channel_names[i]);

	cout << "-------------------------------------------------------------------------------" << endl;


	//Tell the name of the obervables
	RooArgSet *obs = (RooArgSet*) mc->GetObservables();
	obs->Print();

	std::vector<RooRealVar*> x_vector;
	std::vector<RooRealVar*> y_vector;
	std::vector<RooRealVar*> z_vector;

	for(Int_t i=0; i<nchannels;i++)
	{
		x_vector.push_back((RooRealVar*) obs->find(("obs_z_RLc_kinematic_"+channel_names[i]).c_str()));
		y_vector.push_back((RooRealVar*) obs->find(("obs_y_RLc_kinematic_"+channel_names[i]).c_str()));
		z_vector.push_back((RooRealVar*) obs->find(("obs_x_RLc_kinematic_"+channel_names[i]).c_str()));
	}


	int num_channels = meas.GetChannels().size();
	cout<< " num_channels = "<<num_channels<<endl;
	cout<<endl;

	RooCategory *idx = (RooCategory*) obs->find("channelCat");
	idx->Print();
	RooAbsData *data = (RooAbsData*) w->data("obsData");
	data->Print();
	RooPlot *frame = new RooPlot();
	frame = x_vector[0]->frame();
	data->plotOn(frame, RooFit::DrawOption("ZP"));
	TCanvas *c = new TCanvas("c","c");
	frame->Draw();
	c->SaveAs("plot.png");


	mc->GetNuisanceParameters()->Print("v");
	int  nPlotsMax = 1000;
	cout <<" check expectedData by category"<<endl;
	RooDataSet* simData=NULL;
	RooSimultaneous* simPdf = NULL;
	if(strcmp(mc->GetPdf()->ClassName(),"RooSimultaneous")==0){
		cout <<"Is a simultaneous PDF"<<endl;
		simPdf = (RooSimultaneous *)(mc->GetPdf());
	} else {
		cout <<"Is not a simultaneous PDF"<<endl;
	}

	RooAbsPdf *pdf = model_->getPdf(idx->getLabel());
	cout<<"***********************"<<endl;
	cout<<endl;
	cout<<"   INTEGRAL PDF BEFORE FIT: "<< endl;
	//pdf->getNormIntegral(*obs)->getVal()<<endl;
	cout<<pdf->createIntegral(*obs)->getVal()<<endl;
	cout<<endl;
	cout<<"***********************"<<endl;

	RooSimultaneous *model = new RooSimultaneous("simPdf_modified","simPdf_modified",*idx);
	RooAbsPdf* pdf_;
	RooCustomizer* cust;
	//RooRealVar *alpha_mu_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_mu_shape_unc");
	//RooRealVar *alpha_tau_shape_unc = (RooRealVar*) mc->GetNuisanceParameters()->find("alpha_tau_shape_unc");


	//RooRealVar* alpha = new RooRealVar("alpha","alpha", start_parameters["pha"][0], start_parameters["pha"][1], start_parameters["pha"][2]);

	for (Int_t i =0;i< nchannels;i++)
	{
		cout<<i<<endl;
		idx->setIndex(i);
		pdf_=model_->getPdf(idx->getLabel());
		cust = new RooCustomizer(*pdf_,"cust");
		//cust->replaceArg(*alpha_mu_shape_unc,*alpha);
		//cust->replaceArg(*alpha_tau_shape_unc,*alpha);
		model->addPdf((RooAbsPdf&)*cust->build(),idx->getLabel());
	}

	RooStats::HistFactory::HistFactorySimultaneous* model_hf = new RooStats::HistFactory::HistFactorySimultaneous(*model);
	model_hf->indexCat().Print();

	RooAbsReal* nll_hf = model_hf->createNLL(*data,RooFit::Offset(kTRUE));//,RooFit::NumCPU(8));

	RooMinuit* minuit_hf = new RooMinuit(*nll_hf);
	RooArgSet *temp = new RooArgSet();

	nll_hf->getParameters(temp)->Print("V");

	minuit_hf->setErrorLevel(0.5);
	minuit_hf->setStrategy(2);
	minuit_hf->fit("smh");

	RooFitResult *fitResult=minuit_hf->save("TempResult","TempResult");

 for(Int_t i=0; i<nchannels;i++)
        {
		if(channel_names[i]=="Isolated")
			blindResult(fitResult);
	}
	
	std::cout <<"-------CHECK---------------------------------" <<fitResult->edm() << std::endl;
	//Verbose printing: Basic info, values of constant parameters, initial and
    // final values of floating parameters, global correlations
	fitResult->Print("V");
	
	//Summary printing: Basic info plus final values of floating fit parameters
	//fitResult->Print();
for(Int_t i=0; i<nchannels;i++)
        { 
	idx->setIndex(i);
		PlotFrame(x_vector[i],"M_{miss}^{2}",data,model,idx,-2,14,"[GeV^{2}/c^{4}]",channel_names[i],kTRUE);
	PlotFrame(y_vector[i],"E_{l}",data,model,idx,0,2600,"[MeV/c^{2}]",channel_names[i]);
	PlotFrame(z_vector[i],"q^{2}",data,model,idx,-2,14,"[GeV^{2}/c^{4}]",channel_names[i]);
	/*	
	PlotFrame(x_vector[0],"M_{miss}^{2}",data,model_hf,idx,-2,14,"[GeV^{2}/c^{4}]",1);
	PlotFrame(y_vector[0],"E_{l}",data,model_hf,idx,0,2600,"[MeV/c^{2}]");
	PlotFrame(z_vector[0],"q^{2}",data,model_hf,idx,-2,14,"[GeV^{2}/c^{4}]");
	*/
	}

// Access list of final fit parameter values
	if(fitResult->status()==0)
		cout<<" Fit converged "<<endl;
	else
		cout<<"!!!! ATTENTION !!!!    Fit NON CONVERGING!  "<<endl;

    cout << "final value of floating parameters" << endl ;
    fitResult->floatParsFinal().Print("s") ;

	
	//RooRealVar* RLc_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLtau");
//Double_t fitOut[] = {RLc_fitresult->getVal(),RLc_fitresult->errorVar()->getVal()};

//	vector<Double_t> mean_and_error(fitOut,fitOut+sizeof(fitOut)/sizeof(Double_t));
	return fitResult;
}



void SysFit::blindResult(RooFitResult *fitResult)
{
    RooRealVar* RLc_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLctau_Isolated");
    RLc_fitresult->setRange(-999999,999999);
    Double_t true_val = RLc_fitresult->getVal();
    RLc_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha).Uniform(0,100)))*TRandom3(beta).Uniform(0,100)*true_val+TRandom3(gamma).Uniform(0,100));
    RooRealVar* RLcStar_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLcstartau_Isolated");
    RLcStar_fitresult->setRange(-999999,999999);
    RLcStar_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha_s).Uniform(0,100)))*TRandom3(beta_s).Uniform(0,100)*RLcStar_fitresult->getVal()+TRandom3(gamma_s).Uniform(0,100));
}
