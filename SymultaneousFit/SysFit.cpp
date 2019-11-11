#include "SysFit.h"
#include <TStyle.h>
#include <TRandom3.h>

using namespace std;
using namespace RooStats;
using namespace RooFit;
using namespace HistFactory;

SysFit::SysFit()
{
	y_misid_cf = 5.E3;
	y_comb_cf = 7.E3;
	w_misid_cf = 29.;
	w_comb_cf = 54.;
	alpha = 17;
	beta = 34;
	gamma = 98;
	alpha_s = 456;
	beta_s = 172;
	gamma_s = 23;
	BBeast = false;
	ShapeUnc = false;
	GaussConstr = false;
}

Bool_t SysFit::ActivateShapeUncertainties(string name,RooStats::HistFactory::Channel* chan)
{
	if (chan->GetName()=="RLc_kinematic_Kenriched")
	{
		if(name=="2charm")
			ShapeUnc=true;
	}
	if (chan->GetName()=="RLc_kinematic_Isolated")
	{
		if(name=="2charm"||name=="tau"||name=="mu")
			ShapeUnc=true;
	}
	return ShapeUnc;
}

Bool_t SysFit::ActivateGaussianConstraint(string name,RooStats::HistFactory::Channel* chan)
{
	if (chan->GetName()=="RLc_kinematic_Kenriched")
	{
		if(name=="MISID"||name=="Combinatorial")
			GaussConstr=false;
	}
	return GaussConstr;
}

vector <string> SysFit::NameChannels()
{
	vector <string> channel_names;
	//channel_names.push_back(string("Isolated"));
	channel_names.push_back(string("Kenriched"));
	return channel_names;
}

map<string,vector<Double_t>> SysFit::GetStartParameters(string ch_name)
{
	map<string,vector<Double_t>> start_parameters;
	if(ch_name=="Isolated")
	{
		start_parameters.insert(pair<string, vector<Double_t>>("mu",{1.05E6,3.E3,3.E7}));
		start_parameters.insert(pair<string, vector<Double_t>>("2charm",{2.02E5,3.E3,3.E7}));
		start_parameters.insert(pair<string, vector<Double_t>>("starmu",{3.38E5,3.E3,3.E7}));
		start_parameters.insert(pair<string, vector<Double_t>>("MISID",{7.47E4,10,3.E6}));
		start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{5.56E4,1.E1,3.E6}));
		start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.05,1E-9,1.}));
		start_parameters.insert(pair<string, vector<Double_t>>("startau",{0.05,1E-9,1.}));
		start_parameters.insert(pair<string, vector<Double_t>>("pha",{0,-3,3}));
	}
	else if (ch_name=="Kenriched")
	{
		
		start_parameters.insert(pair<string, vector<Double_t>>("mu",{1.4E3,100,1.E4}));
		start_parameters.insert(pair<string, vector<Double_t>>("starmu",{1.E3,10,1.E4}));
		start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.,1E-9,1.}));
		start_parameters.insert(pair<string, vector<Double_t>>("2charm",{1.E3,10,1.E4}));
		start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{3.2E3,1.E3,5.E3}));
		start_parameters.insert(pair<string, vector<Double_t>>("MISID",{4.6E3,100,1.E4}));

		/*
		start_parameters.insert(pair<string, vector<Double_t>>("mu",{1.4E3,1.E3,2.E3}));
		start_parameters.insert(pair<string, vector<Double_t>>("starmu",{1.E3,500,2.E3}));
		start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.,1E-9,1E-4}));
		start_parameters.insert(pair<string, vector<Double_t>>("2charm",{1.E3,500,2.E3}));
		start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{7.E3,5.E3,9.E3}));
		start_parameters.insert(pair<string, vector<Double_t>>("MISID",{5E3,3E3,7.E3}));

		start_parameters.insert(pair<string, vector<Double_t>>("mu",{1E3,100,1.E4}));
                start_parameters.insert(pair<string, vector<Double_t>>("starmu",{4.E3,10,1.E4}));
                start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.,1E-9,1.}));
                start_parameters.insert(pair<string, vector<Double_t>>("2charm",{3.E3,10,1.E4}));
                start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{6.E3,1.E3,1.E4}));
                start_parameters.insert(pair<string, vector<Double_t>>("MISID",{2.E3,100,1.E4}));
		*/
	}
	
	return start_parameters;
}


void SysFit::SetStartParameters(map<string,vector<Double_t> > parameters)
{
	for(std::map<std::string,std::vector<Double_t> > ::iterator iP= parameters.begin(); iP != parameters.end();++iP)
	{
		start_parameters[iP->first.substr(2)][0] = iP->second[0];
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

vector<string> SysFit::GetCategory(map<string,vector<Double_t>> start_parameters)
{       
	vector<string> category;
        // declaring iterators
        map<string, vector<Double_t>>::iterator it1 ;
        for (it1 = start_parameters.begin(); it1!=start_parameters.end(); it1++)
        {       
                category.push_back(string(it1->first));
        }

	return category;
}

TString SysFit::GetComponentName(TString component)
{
	TString name = "";
	if (component.Contains("h_mu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #mu #nu";
	else if(component.Contains("h_starmu_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}^{*} #mu #nu_{#mu}";
	else if(component.Contains("h_startau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c}^{*} #tau #nu_{#tau}";
	else if(component.Contains("h_tau_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} #tau #nu_{#tau}";
	else if(component.Contains("h_2charm_")) name = "#Lambda_{b} #rightarrow #Lambda_{c} X_{c}";
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
	else if(component.Contains("h_tau_"))color = kRed;
	else if(component.Contains("h_2charm_")) color = kGreen;
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
	if(shapeUncert && type.find("2charm_mbody") !=std::string::npos)
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
	if(type=="tau")
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
		if(type=="Combinatorial")
			sample.AddOverallSys("NormUnc_comb_"+category,1.,w_comb_cf);
		if(type=="MISID")
			sample.AddOverallSys("NormUnc_misid_"+category,1.,w_misid_cf);

	}

	if(BBeast) sample.ActivateStatError();
	(*chan)->AddSample(sample);
}

void SysFit::PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,RooStats::HistFactory::HistFactorySimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, bool legend=kFALSE)
{
        kinemObserv->setRange(title,plotStart,plotEnd); //define the variable you want to plot

        //Frame for plotting the fit
        RooPlot *frame = new RooPlot();
        frame = kinemObserv->frame(RooFit::Title(title));
        //Frame for plotting the pulls
        RooPlot *pframe = new RooPlot();
        pframe = kinemObserv->frame("");

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
	while (comp = (TObject*) comp_it->Next())
        {
                component_name = comp->GetName();
                if (comp->InheritsFrom("RooRealSumPdf") && (TString(component_name).Contains("model")))
                {
                        //cout<<component_name<<endl;
                        model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
                }
        }

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
        data->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::DataError( RooAbsData::Poisson), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));

        //cout<<"******************************************"<<endl;
        for (int i = 0; i < active_components_names.size(); ++i)
        {
                std::cout << active_components_names[i] << std::endl;
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
                model_sumpdf->plotOn(frame, RooFit::Slice(*idx),RooFit::ProjWData(*idx,*data),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));
        }
        //cout<<"******************************************"<<endl;

        //--------Plot the sum of the backgrounds
        model_sumpdf->plotOn(frame, RooFit::Slice(*idx),RooFit::ProjWData(*idx,*data),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::FillStyle(3004),RooFit::FillColor(kMagenta-3),RooFit::Components(TotComponents)), RooFit::MoveToBack(),RooFit::Range("title");


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

        pframe->SetMaximum(7.);
        pframe->SetMinimum(-7.);

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

        c1->Print(TString("plots/Fit_")+GetFitVarName(title)+TString(".pdf"));

}



vector<Double_t> SysFit::Fit()
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
		map<string,vector<Double_t>> start_parameters = GetStartParameters(channel_names[i]);
		vector<string> category = GetCategory(start_parameters);
		string filename = string("RootFiles/Histos_")+channel_names[i]+string(".root");

		//Add the samples to the channels
		for(Int_t j=0; j<category.size(); j++)
		{
			cout<<category[j]<<endl;
			if(category[j]=="pha")
				continue;
			ShapeUnc = ActivateShapeUncertainties(category[j],channels[i]);
			GaussConstr = ActivateGaussianConstraint(category[j],channels[i]);
			if(category[j]!="tau"&&category[j]!="2charm")
				AddSample(category[j],filename,0,GaussConstr,BBeast,&channels[i],start_parameters[category[j]]);
			if(category[j]=="2charm")
			{
				AddSample(category[j]+"_2body",filename,0,GaussConstr,BBeast,&channels[i],start_parameters[category[j]]);
				AddSample(category[j]+"_mbody",filename,ShapeUnc,GaussConstr,BBeast,&channels[i],start_parameters[category[j]]);

			}
			if(category[j]=="tau")//NB:: To modify the 0 adding ShapeUnc
				AddSample(category[j],filename,0,GaussConstr,BBeast,&channels[i],start_parameters[category[j]],start_parameters["mu"]);
		}
		channels[i]->SetData("h_data",filename);
		Data ch_data = channels[i]->GetData();
		ch_data.SetName("h_data_"+channel_names[i]);
		cout<<")****** "<<ch_data.GetName()<<endl;
		meas.AddChannel(*channels[i]);
		//Define the parameter of interest
		meas.SetPOI("RLctau_"+channel_names[i]);
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

	//build a workspace with a pdf and a modelconfig
	RooWorkspace *w = RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
	cout<<"--------------------- Workspace created -------------------"<<endl;
	w->Print();
	cout<<"-----------------------------------------------------------"<<endl;

	//Get the model manually
	RooStats::ModelConfig *mc = (RooStats::ModelConfig*) w->obj("ModelConfig");
	RooSimultaneous *model_ = (RooSimultaneous*)mc->GetPdf();

	//Fix the MC normalisation of the histos
	cout << "-------------------------   Fixing parameters   ------------------------------" << endl;
	//mc->GetNuisanceParameters()->Print();
	((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);
	for(Int_t i=0; i<nchannels; i++)
	{
		//cout<<channel_names[i]<<endl;
		if(channel_names[i]=="Kenriched")
		{
			((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel_names[i]).c_str())))->setVal(1E-15);
			((RooRealVar*)(mc->GetParametersOfInterest()->find(("RLctau_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		}
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcmu_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNctau_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstarmu_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		if(channel_names[i]=="Isolated")
		{
			((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcstartau_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		}
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm_mbody_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNc2charm_2body_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcMISID_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		((RooRealVar*)(mc->GetNuisanceParameters()->find(("mcNcCombinatorial_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		if(ActivateGaussianConstraint("MISID",channels[i]))
			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcMISID_"+channel_names[i]).c_str())))->setConstant(kTRUE);
		if(ActivateGaussianConstraint("NcCombinatorial",channels[i]))
			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcCombinatorial_"+channel_names[i]).c_str())))->setConstant(kTRUE);

			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcMISID_"+channel_names[i]).c_str())))->setConstant(kTRUE);

			((RooRealVar*)(mc->GetNuisanceParameters()->find(("NcCombinatorial_"+channel_names[i]).c_str())))->setConstant(kTRUE);
	}

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
	/*
	//Exchange x,z is for historical reasons
	RooRealVar *x = (RooRealVar*) obs->find("obs_z_RLc_kinematic"); //Mmiss2
	RooRealVar *y = (RooRealVar*) obs->find("obs_y_RLc_kinematic"); //E*_mu
	RooRealVar *z = (RooRealVar*) obs->find("obs_x_RLc_kinematic"); //q2
	x->SetTitle("m^{2}_{miss}");
	x->setUnit("GeV^{2}");
	y->SetTitle("E_{#mu}");
	y->setUnit("MeV");
	z->SetTitle("q^{2}");
	z->setUnit("GeV^{2}");
	*/
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
	// blindResult();
	std::cout <<"-------CHECK---------------------------------" <<fitResult->edm() << std::endl;
	fitResult->Print("V");

PlotFrame(x_vector[0],"M_{miss}^{2}",data,model_hf,idx,-2,14,"[GeV^{2}/c^{4}]",1);
        PlotFrame(y_vector[0],"E_{l}",data,model_hf,idx,0,2600,"[MeV/c^{2}]");
        PlotFrame(z_vector[0],"q^{2}",data,model_hf,idx,-2,14,"[GeV^{2}/c^{4}]");


        cout<<"***********************"<<endl;
        cout<<endl;
        cout<<"   INTEGRAL PDF AFTER FIT: "<< endl;
	//pdf->getNormIntegral(*obs)->getVal()<<endl;
	cout<<pdf->createIntegral(*obs)->getVal()<<endl;
        cout<<endl;
        cout<<"***********************"<<endl;

	vector<double_t> r;
	return r;
}


