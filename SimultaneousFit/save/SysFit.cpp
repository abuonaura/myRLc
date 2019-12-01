#include "SysFit.h"
#include <TStyle.h>
#include <TRandom3.h>

using namespace std;
using namespace RooStats;
using namespace RooFit;
using namespace HistFactory;

//namespace RLC
//{
	SysFit::SysFit(){}
	SysFit::~SysFit(){}

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

	void SysFit::PlotFrame(RooRealVar* kinemObserv,const char* title,RooAbsData* data,HistFactorySimultaneous* model, RooCategory* idx,Double_t plotStart, Double_t plotEnd, const char* units, bool legend=kFALSE)
	{
		kinemObserv->setRange(title,plotStart,plotEnd); //define the variable you want to plot

		//Frame for plotting the fit
		RooPlot *frame = new RooPlot();
		frame = kinemObserv->frame(Title(title));
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
				model_sumpdf = (RooRealSumPdf*) components_set->find(component_name);
			}
		}
		//---Take all the components with the name that starts with L_x
		RooArgSet* active_components = model_sumpdf->getComponents();
		vector<TString> active_components_names;

		TIterator* active_comp_it  = active_components->createIterator();
		TObject* active_comp;

		while (active_comp = (TObject*) active_comp_it->Next())
		{
			TString name(active_comp->GetName());
			if(name.Contains("L_x"))
			{
				active_components_names.push_back(name);
			}
		}

		TString TotComponents = "";

		//Plot the data points
		data->plotOn(frame, DrawOption("ZP"), DataError( RooAbsData::Poisson), Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));

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

			//------- Plot each background component
			//model_sumpdf->plotOn(frame, RooFit::Slice(*idx),RooFit::ProjWData(*idx,*data),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(component),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));

			//------- Plot each background component stacking them up
			model_sumpdf->plotOn(frame, RooFit::Slice(*idx),RooFit::ProjWData(*idx,*data),RooFit::ProjectionRange("LbMuNu"),RooFit::DrawOption("F"),RooFit::LineColor(color),RooFit::FillColor(color),RooFit::Components(TotComponents),RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));
		}

		//--------Plot the sum of the backgrounds
		model_sumpdf->plotOn(frame, Slice(*idx),ProjWData(*idx,*data),ProjectionRange("LbMuNu"),DrawOption("F"),FillStyle(3004),FillColor(kMagenta-3),Components(TotComponents)), MoveToBack(),Range("title");


		//-----Construct the pulls plot
		RooHist *pulls = frame->pullHist();
		pulls->SetFillColor(kGray+1);
		pulls->SetLineColor(kWhite);
		pulls->SetMarkerSize(0.01);
		pframe->addPlotable(pulls, "B");


		//Build the legend
		TCanvas *clegend = new TCanvas("clegend", "clegend");
		TLegend leg(0.1,0.1,0.9,0.9);

		if (strcmp(title,"q^{2}")==0)
		{
			vector<TGraph*> component_graphs;
			Int_t CharmComp=0;
			for (int i=0; i<active_components_names.size(); ++i)
			{
				Int_t c = frame->numItems() -3 -i;

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

	map<string,vector<Double_t>> SysFit::GetStartParameters(string ch_name)
	{
		map<string,vector<Double_t>> start_parameters;
		if(strcmp(ch_name,"Isolated")==0)
		{
			start_parameters.insert(pair<string, vector<Double_t>>("mu",{1.05E6,3.E3,3.E7}));
			start_parameters.insert(pair<string, vector<Double_t>>("Ds",{2.02E5,3.E3,3.E7}));
			start_parameters.insert(pair<string, vector<Double_t>>("starmu",{3.38E5,3.E3,3.E7}));
			start_parameters.insert(pair<string, vector<Double_t>>("MISID",{7.47E4,10,3.E6}));
			start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{5.56E4,1.E1,3.E6}));
			start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.05,1E-9,1.}));
			start_parameters.insert(pair<string, vector<Double_t>>("startau",{0.05,1E-9,1.}));
			start_parameters.insert(pair<string, vector<Double_t>>("pha",{0,-3,3}));
		}
		else if (strcmp(ch_name,"Kenriched")==0)
		{

			start_parameters.insert(pair<string, vector<Double_t>>("mu",{1.4E3,100,1.E4}));
			start_parameters.insert(pair<string, vector<Double_t>>("2charm",{1.E3,10,1.E4}));
			start_parameters.insert(pair<string, vector<Double_t>>("tau",{0.,1E-9,1.}));
			start_parameters.insert(pair<string, vector<Double_t>>("starmu",{1.E3,10,1.E4}));
			start_parameters.insert(pair<string, vector<Double_t>>("MISID",{5.E3,100,1.E4}));
			start_parameters.insert(pair<string, vector<Double_t>>("Combinatorial",{7.E3,1.E3,1.E4}));
		}
		return start_parameters;
	}


	void SysFit::AddSample(string type, string inputFile, bool shapeUncert, bool GaussConstraints, const bool BBeast, Channel* chan, vector<Double_t> params, vector<Double_t> mu_params={0})
	{
		//Define the sample
		Sample sample;
		//set the file from where to take the template
		sample.SetInputFile(inputFile);
		//get the normalisation factor to 1
		Double_t mcNorms = GetHistoNormalisation(inputFile, type);

		//Code for shape variations for 2charm sample
		if(shapeUncert && type=="2charm")
		{
			sample.SetName("h_"+type);
			sample.SetHistoName("h_"+type);
			sample.AddHistoSys(type+"_variation","h_"+type+"_1ml",inputFile,"", "h_"+type+"_1pl",inputFile,"");

			sample.AddHistoSys(type+"_quadratic_variation","h_"+type+"_1mq",inputFile,"", "h_"+type+"_1pq",inputFile,"");
		}
		if(shapeUncert && type!="2charm")
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
			sample.AddNormFactor("Ncmu",mu_params[0],mu_params[1],mu_params[2]);
			sample.AddNormFactor("RLc"+type, params[0],params[1],params[2]);
		}
		else
			sample.AddNormFactor("Nc"+type, params[0],params[1],params[2]);
		sample.AddNormFactor("mcNc"+type,mcNorms,1e-20,1);

		//----- Add Gaussian Constraint for MISID/Combinatorial for control fit part
		if (GaussConstraints)
		{
			if(type=="Combinatorial")
				sample.AddOverallSys("NormUnc_comb",1.,w_comb_cf);
			if(type=="MISID")
				sample.AddOverallSys("NormUnc_misid",1.,w_misid_cf);

		}

		if(BBeast) sample.ActivateStatError();
		chan->AddSample(sample);

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




	//---> To add to .h
	vector <string> SysFit::NameChannels()
	{
		vector <string> channel_names;
		channel_names.push_back(string("Isolated"));
		channel_names.push_back(string("Kenriched"));
		return channel_names;
	}

	vector<Double_t> SysFit::Fit()
	{

		//Define the measurement
		Measurement meas("RLc","Rb");
		//Define where to save the workspace
		meas.SetOutputFilePrefix("results/RLc");
		//Tell histfactory not to run the fit itself
		meas.SetExportOnly(kTRUE);
		//Define the parameter of interest
		meas.SetPOI("RLctau");
		//Set the luminosity
		meas.SetLumi(1.0);
		meas.SetLumiRelErr(0.05);

		//?
		RooRealVar* alpha = new RooRealVar("alpha","alpha", start_parameters["pha"][0], start_parameters["pha"][1], start_parameters["pha"][2]);

		//Define a channel for each category
		vector <string> channel_names = NameChannels();
		Int_t nchannels = channel_names.size();
		vector <Channel*> channels;
		for(Int_t i=0; i<nchannels; i++)
		{
			channels.push_back(new Channel((string("RLc_kinematic_")+channel_names[i]).c_str()));
			channels[i]->SetStatErrorConfig(1e-5, "Poisson");
			map<string,vector<Double_t>> start_parameters = GetStartParameters(channel_names[i]);
			cout<<endl;
			cout<<"---------------------------------------------------------"<<endl;
			PrintStartParams(channel_names[i],start_parameters);
			cout<<"---------------------------------------------------------"<<endl;
			cout<<endl;
		}





	}










	//------------------> Not modified!
	std::map<std::string,std::vector<Double_t> > SysFit::GetFitResult()
	{
		RooArgList floatParams = (RooArgList)fitResult->floatParsFinal();
		for(Int_t iName=0;iName<floatParams.getSize();iName++)
		{
			RooRealVar* fit_parameter = (RooRealVar*)floatParams.find(floatParams[iName].GetTitle());
			std::string name(fit_parameter->GetTitle());
			fit_result[(std::string)fit_parameter->GetTitle()] = {fit_parameter->getVal()};
		}
		return fit_result;
	}

	void SysFit::SetStartParameters(std::map<std::string,std::vector<Double_t> > parameters)
	{
		for(std::map<std::string,std::vector<Double_t> > ::iterator iP= parameters.begin(); iP != parameters.end();++iP)
		{
			start_parameters[iP->first.substr(2)][0] = iP->second[0];
		}
	}

	void SysFit::RLcPlot()
	{ 
		Int_t n = 2;
		Double_t RLc_value[n],RLcStar_value[n],RLc_error[n],RLcStar_error[n], angle[n];
		this->Fit();
		RooRealVar* RLctau = (RooRealVar*)fitResult->floatParsFinal().find("RLtau");
		RooRealVar* RLcLcStartau = (RooRealVar*)fitResult->floatParsFinal().find("RLLcStartau");
		angle[0] = correlation;
		RLc_value[0] = RLctau->getVal();
		RLc_error[0] = RLctau->errorVar()->getVal();
		RLcStar_value[0] = RLcLcStartau->getVal();
		RLcStar_error[0] = RLcLcStartau->errorVar()->getVal();
		std::map<std::string,std::vector<Double_t> > fitRes = this->GetFitResult();
		this->SetStartParameters(fitRes);
		this->SetBB(1);
		this->Fit();
		RLctau = (RooRealVar*)fitResult->floatParsFinal().find("RLtau");
		RLcLcStartau = (RooRealVar*)fitResult->floatParsFinal().find("RLLcStartau");
		angle[1] = correlation;
		RLc_value[1] = RLctau->getVal();
		RLc_error[1] = RLctau->errorVar()->getVal();
		RLcStar_value[1] = RLcLcStartau->getVal();
		RLcStar_error[1] = RLcLcStartau->errorVar()->getVal();
		gStyle->SetCanvasPreferGL(1);
		gStyle->SetOptStat(0);
		TCanvas *cRLc = new TCanvas("RLc","RLc",200,10,500,300);
		// cRLc->Range(0,0,0.15,0.15);

		TH2F *frame = new TH2F("frame", "error ellipse",
				1, 0, 0.1, 
				1, 0, 0.15 );
		frame->GetXaxis()->SetTitle("R(L_{c})");
		frame->GetYaxis()->SetTitle("R(L_{c}*)");
		frame->Draw();

		TH1F *toLegendBB = new TH1F("toLegBB","toLegBB",1,-2,2);
		toLegendBB->SetFillColor(7);
		toLegendBB->Draw("same");

		TH1F *toLegend = new TH1F("toLeg","toLeg",1,-2,2);
		toLegend->SetFillColor(5);
		toLegend->Draw("same");

		TLegend *leg = new TLegend(0.7,0.8,.99,.99);

		leg->AddEntry("toLeg","Normal fitting","f");
		leg->AddEntry("toLegBB","Barlow-Beaston","f");
		leg->Draw();

		drawEllipse(RLc_value[1],RLc_error[1],RLcStar_value[1],RLcStar_error[1],7,-45*angle[1]);

		drawEllipse(RLc_value[0],RLc_error[0],RLcStar_value[0],RLcStar_error[0],5,-45*angle[0]);


		cRLc->SaveAs("ellipse.pdf");
		cout<<"Normal fitting: \nMean:    "<<RLc_value[0]<<"      Sigma: "<<RLc_error[0]<<"      correlation:    "<<angle[0]<<endl;
		cout<<"BB fitting: \nMean:    "<<RLc_value[1]<<"      Sigma: "<<RLc_error[1]<<"      correlation:    "<<angle[1]<<endl;
	}

	void SysFit::drawEllipse(Double_t x,Double_t x_err,Double_t y,Double_t y_err,Int_t color, Double_t correlation)
	{
		TEllipse *err_ellipse = new TEllipse(x,y,x_err,y_err,0,360,Int_t(correlation));
		err_ellipse->SetLineWidth(1);
		err_ellipse->SetLineColor(color);
		err_ellipse->SetFillColorAlpha(color,0.3);
		err_ellipse->Draw();

	}

	void SysFit::blindResult()
	{
		RooRealVar* RLc_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLtau");
		RLc_fitresult->setRange(-999999,999999);
		Double_t true_val = RLc_fitresult->getVal();
		RLc_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha).Uniform(0,100)))*TRandom3(beta).Uniform(0,100)*true_val+TRandom3(gamma).Uniform(0,100));
		RooRealVar* RLcStar_fitresult = (RooRealVar*)fitResult->floatParsFinal().find("RLLcStartau");
		RLcStar_fitresult->setRange(-999999,999999);
		RLcStar_fitresult->setVal(std::pow(-1,Int_t(TRandom3(alpha_s).Uniform(0,100)))*TRandom3(beta_s).Uniform(0,100)*RLcStar_fitresult->getVal()+TRandom3(gamma_s).Uniform(0,100));
	}







//}
