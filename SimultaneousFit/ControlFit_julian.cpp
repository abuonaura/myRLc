#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLegend.h"
#include "RooPolynomial.h"
#include "RooCustomizer.h"
#include "RooChi2Var.h"
#include "RooAbsData.h"
#include "RooRealSumPdf.h"
#include "RooPoisson.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
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
#include "TString.h"
#include "TH1D.h"
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
#include <boost/program_options/option.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <fstream>
#include "FitConfig.hh"
#include "FitParams.hh"
#include "RDplus_FitParams.hh"
#include "RDplus_FitParams_Common.hh"
#include "BGLFFParamConstr.hh"
#include "RooTFnBinding.h" 
#include "MCStatPlot.hh"
#include "RooUniform.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "RooHammerModel.h"
#include "Hammer/Hammer.hh"
#include <ctime>

//Generated files that contain the effects of the systematics

TString pretty_component_name(TString component_name)
{
  TString pretty_component_name(""); 
  if(component_name.Contains("L_x_hNMisID_")){pretty_component_name = "MisID";}
  else if(component_name.Contains("L_x_hNWS_")){pretty_component_name = "Combinatorial";}
  else if(component_name.Contains("L_x_hNDpMuNu_")){pretty_component_name = "B^{0}#rightarrow D^{-} #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNDstMuNu_")){pretty_component_name = "B^{0}#rightarrow D*^{-} #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2DststmunuHigher_")){pretty_component_name = "B^{0} #rightarrow (D_{J}** #rightarrow D^{-} X ) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D2star")){pretty_component_name = "B^{0} #rightarrow (D_{2}* #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D0star")){pretty_component_name = "B^{0} #rightarrow (D_{0}* #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D1prime")){pretty_component_name = "B^{0} #rightarrow (D_{1}' #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D1")){pretty_component_name = "B^{0} #rightarrow (D_{1} #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D2star")){pretty_component_name = "B^{#pm} #rightarrow (D_{2}* #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D0star")){pretty_component_name = "B^{#pm} #rightarrow (D_{0}* #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D1prime")){pretty_component_name = "B^{#pm} #rightarrow (D_{1}' #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D1")){pretty_component_name = "B^{#pm} #rightarrow (D_{1} #rightarrow D^{-} X) #mu^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNDpTauNu_")){pretty_component_name = "B^{0} #rightarrow D^{-} #tau^{+}#nu_{#mu}";}
  else if(component_name.Contains("L_x_hNDstTauNu_")){pretty_component_name = "B^{0} #rightarrow D*^{-} #tau^{+} #nu_{#mu}";}
  else if(component_name.Contains("L_x_hNBd2DD_DDs_")){pretty_component_name = "B^{0} #rightarrow D^{-} (Xc #rightarrow #mu^{+} #nu_{#mu} X)";}
  else if(component_name.Contains("L_x_hNBu2DD_DDs_")){pretty_component_name = "B^{#pm} #rightarrow D^{-} (Xc #rightarrow #mu^{+} #nu_{#mu} X)";}
  else if(component_name.Contains("L_x_hNBdDDs_")){pretty_component_name = "B^{0} #rightarrow D^{-} (Xc #rightarrow #tau^{+} #nu_{#tau} X)";}
  else if(component_name.Contains("L_x_hNBuDDs_")){pretty_component_name = "B^{#pm} #rightarrow D^{-} (Xc #rightarrow #tau^{+} #nu_{#tau} X)";}
  else if(component_name.Contains("L_x_hNBd2DD_MultiBody_")){pretty_component_name = "B^{0}#rightarrow D^{-} (Xc #rightarrow #mu^{+} #nu_{#mu} X) X";}
  else if(component_name.Contains("L_x_hNBu2DD_MultiBody_")){pretty_component_name = "B^{#pm} #rightarrow D^{-} (Xc #rightarrow #mu^{+} #nu_{#mu} X) X";}
  else if(component_name.Contains("L_x_hNLbDLc_")){pretty_component_name = "#Lambda_{b} #rightarrow D^{-} (#Lambda_{c} #rightarrow #mu^{+} #nu_{#mu} X ) X";}
  else{pretty_component_name = component_name;}
  
  pretty_component_name = "#font[12]{"+pretty_component_name+"}";

  return pretty_component_name; 
}

int component_color(TString component_name)
{
  int color; 
  if(component_name.Contains("L_x_hNMisID_")){color = kGray;}
  else if(component_name.Contains("L_x_hNWS_")){color = kGreen+2;}
  else if(component_name.Contains("L_x_hNDpMuNu_")){color = 4;}
  else if(component_name.Contains("L_x_hNDstMuNu_")){color = 38;}
  else if(component_name.Contains("L_x_hNBd2DststmunuHigher_")){color = kViolet+1;}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D2star")){color = kMagenta+2;}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D0star")){color = kMagenta+2;}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D1")){color = kMagenta+2;}
  else if(component_name.Contains("L_x_hNBd2Dststmunu_D1prime")){color = kMagenta+2;}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D2star")){color = kMagenta;}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D0star")){color = kMagenta;}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D1")){color = kMagenta;}
  else if(component_name.Contains("L_x_hNBu2Dststmunu_D1prime")){color = kMagenta;}
  else if(component_name.Contains("L_x_hNDpTauNu_")){color = 2;}
  else if(component_name.Contains("L_x_hNDstTauNu_")){color=46;}
  else if(component_name.Contains("L_x_hNBd2DD_DDs_")){color = kYellow;}
  else if(component_name.Contains("L_x_hNBu2DD_DDs_")){color = kYellow-3;}
  else if(component_name.Contains("L_x_hNBdDDs_")){color = kRed-5;}
  else if(component_name.Contains("L_x_hNBuDDs_")){color = kRed-8;}
  else if(component_name.Contains("L_x_hNBd2DD_MultiBody_")){color = kOrange;}
  else if(component_name.Contains("L_x_hNBu2DD_MultiBody_")){color = kOrange-3;}
  else if(component_name.Contains("L_x_hNLbDLc_")){color = kCyan-3;}
  else{color = kBlack;} //Draw in back the components you have not found
  
  return color; 
}

namespace po = boost::program_options;

void LoadStartValues(RooStats::ModelConfig* mc, RDplus_common_parameters* common_parameters, std::vector<RDplus_parameters*>& toyparams_vector){
  
  parameter* param;

  for (unsigned int k = 0 ; k<common_parameters->nparameters(); ++k){
    param = common_parameters->get_parameter(k);

    RooRealVar* model_param_nuisance  = (RooRealVar*) mc->GetNuisanceParameters()->find((param->get_name()).c_str());
    RooRealVar* model_param_poi       = (RooRealVar*) mc->GetParametersOfInterest()->find((param->get_name()).c_str());
    if (model_param_nuisance){
        model_param_nuisance->setVal(param->get_value());
    }
    else if(model_param_poi){
        model_param_poi->setVal(param->get_value());
    }
    else{
      std::cout << "INFO: Parameter "<< param->get_name() << " not found in the model. Leaving as it is while loading the start values!" << std::endl;
    }
    //If the parameter has not been found in the model, then leave it as it is
  }

  for (unsigned int i=0; i<toyparams_vector.size(); i++){
    for (unsigned int k=0; k<toyparams_vector[i]->nparameters(); ++k){
      param = toyparams_vector[i]->get_parameter(k);

      RooRealVar* model_param_nuisance  = (RooRealVar*) mc->GetNuisanceParameters()->find((param->get_name()).c_str());
      RooRealVar* model_param_poi       = (RooRealVar*) mc->GetParametersOfInterest()->find((param->get_name()).c_str());
      if (model_param_nuisance){
        model_param_nuisance->setVal(param->get_value());
      }
      else if(model_param_poi){
        model_param_poi->setVal(param->get_value());
      }
      else{
        std::cout << "INFO: Parameter "<< param->get_name() << " not found in the model. Leaving as it is while loading the start values!" << std::endl;
      }
    }
  }

  return;
}

void saveResult(RooStats::ModelConfig* mc, RDplus_common_parameters* common_parameters, std::vector<RDplus_parameters*>& toyparams_vector){

  parameter* param;

  for (unsigned int k = 0 ; k<common_parameters->nparameters(); ++k){
    param = common_parameters->get_parameter(k);

    RooRealVar* model_param_nuisance  = (RooRealVar*) mc->GetNuisanceParameters()->find((param->get_name()).c_str());
    RooRealVar* model_param_poi       = (RooRealVar*) mc->GetParametersOfInterest()->find((param->get_name()).c_str());
    if (model_param_nuisance){
      param->set_value(model_param_nuisance->getVal());
      param->set_error_up(model_param_nuisance->getError());
      param->set_error_down(model_param_nuisance->getError());
    }
    else if(model_param_poi){
      param->set_value(model_param_poi->getVal());
      param->set_error_up(model_param_poi->getError());
      param->set_error_down(model_param_poi->getError());
    }
    else{
      std::cout << "INFO: Parameter "<< param->get_name() << " not found in the model. Leaving as it is while saving the text files!" << std::endl;
    }
    //If the parameter has not been found in the model, then leave it as it is
  }

  for (unsigned int i=0; i<toyparams_vector.size(); i++){
    for (unsigned int k=0; k<toyparams_vector[i]->nparameters(); ++k){
      param = toyparams_vector[i]->get_parameter(k);

      RooRealVar* model_param_nuisance  = (RooRealVar*) mc->GetNuisanceParameters()->find((param->get_name()).c_str());
      RooRealVar* model_param_poi       = (RooRealVar*) mc->GetParametersOfInterest()->find((param->get_name()).c_str());
      if (model_param_nuisance){
        param->set_value(model_param_nuisance->getVal());
        param->set_error_up(model_param_nuisance->getError());
        param->set_error_down(model_param_nuisance->getError());
      }
      else if(model_param_poi){
        param->set_value(model_param_poi->getVal());
        param->set_error_up(model_param_poi->getError());
        param->set_error_down(model_param_poi->getError());
      }
      else{
        std::cout << "INFO: Parameter "<< param->get_name() << " not found in the model. Leaving as it is while saving the text files!" << std::endl;
      }
    }
  }

  return;
}

void saveFittedYields(char* outputfolder, char* outputfileprefix, std::vector<std::string> name_suffix_vector, std::vector<std::string>& isocats ,RooStats::ModelConfig* mc, RDplus_common_parameters* common_parameters, std::vector<RDplus_parameters*>& toyparams_vector){

  std::string outputfile_name;
  std::string outputfile_common_name;
  std::ofstream fityieldsfile_common;
  std::ofstream fityieldsfile;


  outputfile_common_name = std::string(outputfolder)+std::string("/")+std::string(outputfileprefix)+std::string("_common.txt");
  fityieldsfile_common.open(outputfile_common_name.c_str());


  // Print the parameters in the file
  for (unsigned int k=0; k < common_parameters->nparameters(); k++){
    parameter* param = common_parameters->get_parameter(k);
    fityieldsfile_common <<std::fixed<< param->get_name() << " " <<  param->get_value() << " " << param->get_min() << " " << param->get_max() << std::endl;
  }

  fityieldsfile_common.close();

  for (unsigned int i=0; i<name_suffix_vector.size(); ++i){

    outputfile_name = std::string(outputfolder)+std::string("/")+std::string(outputfileprefix)+std::string("_cat")+std::string(isocats[i])+std::string(".txt");
    fityieldsfile.open(outputfile_name.c_str());

    for (unsigned int k=0; k < toyparams_vector[i]->nparameters(); k++){
      parameter* param = toyparams_vector[i]->get_parameter(k);

      std::string original_name = std::string(param->get_name());
      std::string name = std::string(original_name).erase(original_name.find(name_suffix_vector[i]), original_name.size());

      fityieldsfile <<std::fixed<< name << " " <<  param->get_value() << " " << param->get_min() << " " << param->get_max() << std::endl;
    }

    fityieldsfile.close();
  }

  return;
}

void plotProjection(RooCategory* idx,std::map<std::string,RooDataHist*> data,RooSimultaneous* model,RooRealVar* obs,TString folder, TString fitvar, TString label, RooStats::ModelConfig* mc,RooWorkspace* w, TString output_file_name, TFile* fin, const char* name_suffix, std::map<std::string, HistogramErrors> error_histos, std::string selection="")
{
  
  gROOT->ProcessLine(".x mylhcbStyle.C++");

  //---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
  //RooAbsPdf* model_pdf = new RooProdPdf(("model_pdf_chan"+std::to_string(idx->getIndex())).c_str(),("model_pdf_chan"+std::to_string(idx->getIndex())).c_str(),*model->getPdf(idx->getLabel()),sweightcorr_inv_pdfs.at(idx->getIndex())); // This gets the pdf in the channel defined by the idx.
  RooAbsPdf* model_pdf = model->getPdf(idx->getLabel());
  
  RooArgSet* components_set = model_pdf->getComponents(); 
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
  //RooProdPdf* model_sumpdf = new RooProdPdf(("model_pdf_chan"+std::string(name_suffix)).c_str(),("model_pdf_chan"+std::string(name_suffix)).c_str(),*model_sumpdf_,sweightcorr_inv_pdfs.at(idx->getIndex()));

  RooPlot *frame = new RooPlot();
  frame->SetName(TString(obs->GetName())+TString("_frame_")+idx->getLabel());
  frame = obs->frame("");
  RooPlot *pframe = new RooPlot();
  pframe->SetName(TString(obs->GetName())+TString("_pframe_")+idx->getLabel());
  pframe = obs->frame(""); //frame used to plot the pulls

  frame->GetXaxis()->SetTitle(label);
  pframe->GetXaxis()->SetTitle(label);

  //---Take all the components with the name that starts with L_x
  RooArgSet* active_components = model_sumpdf->getComponents();

  std::vector<TString> active_components_names;

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
 
  //------Set all the Barlow Beeston parameters back to 1 to see the actual fit result and not to fool yourself! (See Phoebe's Tutorial and Twiki Page)
  RooRealVar *tempVar; //dummy pointer
  TIterator *itergamma = mc->GetNuisanceParameters()->createIterator();
  while((tempVar = (RooRealVar*) itergamma->Next())) {
    TString theName = tempVar->GetName();
    if(theName.Contains("gamma_stat")){
      tempVar->setVal(1.0);
    }
  } 
  //Try to see all the components
  TString components = "";
  //RooAbsData* channelData = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  RooAbsData* channelData = data[idx->getLabel()];
  channelData->plotOn(frame, RooFit::DrawOption("ZP"));
  int colour_count=0;
 
  for (int i = 0; i < active_components_names.size(); ++i)
  {
    std::cout << active_components_names[i] << std::endl;
    if (components != "")
    {
      components+=",";
    }
    components+="*";
    components+=active_components_names[i];
    components+="*";

    int colour_code = component_color(active_components_names[i]);

    //model_sumpdf->plotOn(frame,RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*channelData),RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));
    model_sumpdf->plotOn(frame,RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot").Data()));  
    colour_count++;
  }

  //-----Construct the pulls plot
  //model_sumpdf->plotOn(frame, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*data),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  model_sumpdf->plotOn(frame, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  RooHist *pulls = frame->pullHist();
  pulls->SetFillColor(kGray+1);
  pulls->SetLineColor(kWhite);
  pulls->SetMarkerSize(0.01);
  pframe->addPlotable(pulls, "B");

  //----Construct the MCstat plot
  std::cout << "-----Creating the MCStatPlot------" << std::endl;   
  MCStatPlot* MCplot = new MCStatPlot("MCStatPlot", model_sumpdf);

  MCplot->SetMap(error_histos);
  MCplot->SetHistoNameSuffix(TString(name_suffix));
  MCplot->SetNormalizationMode("data");

  TH1F* MCplot_histogram_projected;

  if(TString(obs->GetName()).Contains("obs_x"))
    MCplot_histogram_projected = (TH1F*) MCplot->GetHistogram1D("x");
  else if(TString(obs->GetName()).Contains("obs_y"))
    MCplot_histogram_projected = (TH1F*) MCplot->GetHistogram1D("y");
  else if (TString(obs->GetName()).Contains("obs_z"))
    MCplot_histogram_projected = (TH1F*) MCplot->GetHistogram1D("z");

  MCplot_histogram_projected->SetFillStyle(3244);
  MCplot_histogram_projected->SetFillColor(9);
  MCplot_histogram_projected->SetMarkerColor(0);
  MCplot_histogram_projected->SetMarkerSize(0.00001);
  std::cout << "----------------------------------" << std::endl;
  
  //data->plotOn(frame, RooFit::DrawOption("ZP"), RooFit::Cut(TString("channelCat==channelCat::")+TString(idx->getLabel())));
  data[idx->getLabel()]->plotOn(frame, RooFit::DrawOption("ZP"));
  TLegend leg(0.1,0.1,0.9,0.9);

  if ((fitvar.Contains("Mmiss2")) || (fitvar.Contains("BDT")))
  {
      std::vector<TGraph*> component_graphs;
      TGraph* graph;
      for (int i=0; i<active_components_names.size(); ++i)
      {
        leg.AddEntry(active_components_names[i]+TString("_plot").Data(), pretty_component_name(active_components_names[i]), "f");  
      }
  }

  //-----------Draw the fit result with the pulls plot below
  TCanvas *clegend = new TCanvas("clegend", "clegend");

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

  Double_t xmin = -5.;
  Double_t xmax = 5.;
  TF1* upperLine_function = new TF1("upperLine", "3.", xmin, xmax);
  TF1* mediumLine_function= new TF1("mediumLine", "0.", xmin, xmax);
  TF1* lowerLine_function = new TF1("lowerLine", "-3.",xmin, xmax);
  RooAbsReal* upperLine_roofit = RooFit::bindFunction(upperLine_function, *obs);
  RooAbsReal* mediumLine_roofit= RooFit::bindFunction(mediumLine_function,*obs);
  RooAbsReal* lowerLine_roofit = RooFit::bindFunction(lowerLine_function, *obs);

  upperLine_roofit->plotOn(pframe, RooFit::LineStyle(3), RooFit::LineColor(12));
  mediumLine_roofit->plotOn(pframe, RooFit::LineStyle(3),RooFit::LineColor(12));
  lowerLine_roofit->plotOn(pframe, RooFit::LineStyle(3), RooFit::LineColor(12));

  pframe->Draw();
  //total_templates_uncertainty.Draw("e2,same");
  MCplot_histogram_projected->Draw("e2,same");

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

  //if ((fitvar.Contains("Mmiss2")) || (fitvar.Contains("BDT"))){
  clegend->cd();
  leg.Draw();
  clegend->Print(folder+TString("/")+TString("Legend.pdf"));
  clegend->Print(folder+TString("/")+TString("Legend.C"));
  //}

  c1->Print(folder+TString("/")+fitvar+TString("_")+output_file_name+TString(name_suffix)+TString(".pdf"));
  c1->Print(folder+TString("/")+fitvar+TString("_")+output_file_name+TString(name_suffix)+TString(".root"));
  c1->Print(folder+TString("/")+fitvar+TString("_")+output_file_name+TString(name_suffix)+TString(".C"));
}

bool IsSampleVetoed(std::string isocat_string, std::string sample_name){

  bool isVetoed = false;

  //---Map that contains the vetoed samples in the corresponding channels
  std::map<std::string,std::vector<std::string>> vetoed_samples;
  
  std::vector<std::string> vector_vetoed_DplusMuNu       {}; //{std::string("10"), std::string("11")};
  std::vector<std::string> vector_vetoed_DstMuNu         {}; //{std::string("10"), std::string("11")};
  std::vector<std::string> vector_vetoed_DplusTauNu      {}; //{std::string("10"), std::string("11"), std::string("2"), std::string("4"), std::string("12"), std::string("13"), std::string("14")};
  std::vector<std::string> vector_vetoed_DstTauNu        {}; //{std::string("10"), std::string("11"), std::string("2"), std::string("4"), std::string("12"), std::string("13"), std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2DD_DDs       {}; //{std::string("10"), std::string("11"), std::string("2") };
  std::vector<std::string> vector_vetoed_Bu2DD_DDs       {}; //{std::string("10"), std::string("11")                   };
  std::vector<std::string> vector_vetoed_BdDDs           {}; //{std::string("10"), std::string("11"), std::string("2"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_BuDDs           {}; //{std::string("10"), std::string("11"), std::string("2"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2DD_MultiBody {}; //{std::string("10"), std::string("11"), std::string("2") };
  std::vector<std::string> vector_vetoed_Bu2DD_MultiBody {}; //{std::string("10"), std::string("11")                   };
  std::vector<std::string> vector_vetoed_Bu2Dstst_D2star {}; //{std::string("10"), std::string("11"), std::string("12")};
  std::vector<std::string> vector_vetoed_Bu2Dstst_D0star {}; //{std::string("10"), std::string("11"), std::string("12")};
  std::vector<std::string> vector_vetoed_Bu2Dstst_D1     {}; //{std::string("10"), std::string("11"), std::string("12")};
  std::vector<std::string> vector_vetoed_Bu2Dstst_D1prime{}; //{std::string("10"), std::string("11"), std::string("12")};
  std::vector<std::string> vector_vetoed_Bd2Dstst_D2star {}; //{std::string("10"), std::string("11"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2Dstst_D0star {}; //{std::string("10"), std::string("11"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2Dstst_D1     {}; //{std::string("10"), std::string("11"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2Dstst_D1prime{}; //{std::string("10"), std::string("11"), std::string("13"),std::string("14")};
  std::vector<std::string> vector_vetoed_Bd2DststHigher  {}; //{std::string("10"), std::string("11")};
  std::vector<std::string> vector_vetoed_LbDLc           {}; //{std::string("10"), std::string("11"),std::string("4"),std::string("13"), std::string("14")};
  std::vector<std::string> vector_vetoed_MisID           {};
  std::vector<std::string> vector_vetoed_WS              {};
  //std::vector<std::string> vector_vetoed_WS              {std::string("4")}; //TODO: Remove this constraint once the new MisID weights are ready

  vetoed_samples[std::string("DplusMuNu")]       = vector_vetoed_DplusMuNu;
  vetoed_samples[std::string("DstMuNu")]         = vector_vetoed_DstMuNu;
  vetoed_samples[std::string("DplusTauNu")]      = vector_vetoed_DplusTauNu;
  vetoed_samples[std::string("DstTauNu")]        = vector_vetoed_DstTauNu;
  vetoed_samples[std::string("Bd2DD_DDs")]       = vector_vetoed_Bd2DD_DDs;
  vetoed_samples[std::string("Bu2DD_DDs")]       = vector_vetoed_Bu2DD_DDs;
  vetoed_samples[std::string("BdDDs")]           = vector_vetoed_BdDDs;
  vetoed_samples[std::string("BuDDs")]           = vector_vetoed_BuDDs;
  vetoed_samples[std::string("Bd2DD_MultiBody")] = vector_vetoed_Bd2DD_MultiBody;
  vetoed_samples[std::string("Bu2DD_MultiBody")] = vector_vetoed_Bu2DD_MultiBody;
  vetoed_samples[std::string("Bu2Dstst_D2star")] = vector_vetoed_Bu2Dstst_D2star;
  vetoed_samples[std::string("Bu2Dstst_D0star")] = vector_vetoed_Bu2Dstst_D0star;
  vetoed_samples[std::string("Bu2Dstst_D1")]     = vector_vetoed_Bu2Dstst_D1;
  vetoed_samples[std::string("Bu2Dstst_D1prime")]= vector_vetoed_Bu2Dstst_D1prime;
  vetoed_samples[std::string("Bd2Dstst_D2star")] = vector_vetoed_Bd2Dstst_D2star;
  vetoed_samples[std::string("Bd2Dstst_D0star")] = vector_vetoed_Bd2Dstst_D0star;
  vetoed_samples[std::string("Bd2Dstst_D1")]     = vector_vetoed_Bd2Dstst_D1;
  vetoed_samples[std::string("Bd2Dstst_D1prime")]= vector_vetoed_Bd2Dstst_D1prime;
  vetoed_samples[std::string("Bd2DststHigher")]  = vector_vetoed_Bd2DststHigher;
  vetoed_samples[std::string("LbDLc")]           = vector_vetoed_LbDLc;
  vetoed_samples[std::string("MisID")]           = vector_vetoed_MisID;
  vetoed_samples[std::string("WS")]              = vector_vetoed_WS;

  std::map<std::string,std::vector<std::string>>::iterator map_it;
  std::vector<std::string>::iterator vec_it;

  map_it = vetoed_samples.find(sample_name);
  if (map_it != vetoed_samples.end()){
    //If you have found the sample in one of the vetoes, check if it is vetoed on this exact channel
    for(vec_it = map_it->second.begin(); vec_it!=map_it->second.end(); ++vec_it){
      if (*vec_it == isocat_string) isVetoed = true; 
    }
  }

  return isVetoed;
}

void plotProjectionInBins(RooCategory* idx,std::map<std::string,RooDataHist*> data,RooSimultaneous* model, RooRealVar* obs,TString folder, TString fitvar, TString label, RooStats::ModelConfig* mc,RooWorkspace* w, TString output_file_name, TFile* fin, const char* name_suffix, std::map<std::string, HistogramErrors> error_histos, std::string selection="")
{
  std::cout << "Inside the function" << std::endl;
  
  gROOT->ProcessLine(".x mylhcbStyle.C++");
  
  //---Get the RooSumPdf. This is needed in order to plot the model from 6.14 onwards
  //RooAbsPdf* model_pdf = new RooProdPdf(("model_pdf_chan"+std::to_string(idx->getIndex())).c_str(),("model_pdf_chan"+std::to_string(idx->getIndex())).c_str(),*model->getPdf(idx->getLabel()),sweightcorr_inv_pdfs.at(idx->getIndex())); // This gets the pdf in the channel defined by the idx.
  RooAbsPdf* model_pdf = model->getPdf(idx->getLabel());
  
  RooArgSet* components_set = model_pdf->getComponents(); 
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
  //RooProdPdf* model_sumpdf = new RooProdPdf(("model_pdf_chan"+std::string(name_suffix)).c_str(),("model_pdf_chan"+std::string(name_suffix)).c_str(),*model_sumpdf_,sweightcorr_inv_pdfs.at(idx->getIndex()));

  RooPlot *frame1 = new RooPlot();
  RooPlot *frame2 = new RooPlot();
  RooPlot *frame3 = new RooPlot();
  RooPlot *frame4 = new RooPlot();

  
  frame1->SetName(TString(obs->GetName())+TString("_frame1_")+idx->getLabel());
  frame2->SetName(TString(obs->GetName())+TString("_frame2_")+idx->getLabel());
  frame3->SetName(TString(obs->GetName())+TString("_frame3_")+idx->getLabel());
  frame4->SetName(TString(obs->GetName())+TString("_frame4_")+idx->getLabel());
  
  frame1 = obs->frame("");
  frame2 = obs->frame("");
  frame3 = obs->frame("");
  frame4 = obs->frame("");

  RooPlot *pframe1 = new RooPlot();
  RooPlot *pframe2 = new RooPlot();
  RooPlot *pframe3 = new RooPlot();
  RooPlot *pframe4 = new RooPlot();

  pframe1->SetName(TString(obs->GetName())+TString("_pframe1_")+idx->getLabel());
  pframe2->SetName(TString(obs->GetName())+TString("_pframe2_")+idx->getLabel());
  pframe3->SetName(TString(obs->GetName())+TString("_pframe3_")+idx->getLabel());
  pframe4->SetName(TString(obs->GetName())+TString("_pframe4_")+idx->getLabel());
  
  pframe1 = obs->frame(""); //frame used to plot the pulls
  pframe2 = obs->frame(""); //frame used to plot the pulls
  pframe3 = obs->frame(""); //frame used to plot the pulls
  pframe4 = obs->frame(""); //frame used to plot the pulls
  
  frame1->GetXaxis()->SetTitle(label);
  frame2->GetXaxis()->SetTitle(label);
  frame3->GetXaxis()->SetTitle(label);
  frame4->GetXaxis()->SetTitle(label);
  pframe1->GetXaxis()->SetTitle(label);
  pframe2->GetXaxis()->SetTitle(label);
  pframe3->GetXaxis()->SetTitle(label);
  pframe4->GetXaxis()->SetTitle(label);

  
  //---Take all the components with the name that starts with L_x
  RooArgSet* active_components = model_sumpdf->getComponents();
  std::vector<TString> active_components_names;

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
  
  //------Set all the Barlow Beeston parameters back to 1 to see the actual fit result and not to fool yourself! (See Phoebe's Tutorial and Twiki Page)
  RooRealVar *tempVar; //dummy pointer
  TIterator *itergamma = mc->GetNuisanceParameters()->createIterator();
  while((tempVar = (RooRealVar*) itergamma->Next())) {
    TString theName = tempVar->GetName();
    if(theName.Contains("gamma_stat")){
      tempVar->setVal(1.0);
    }
  } 
  
  // Get the variable the projection is defined on
  RooRealVar* q2 = (RooRealVar*) mc->GetObservables()->find((std::string("obs_x_RDplus_kinematic")+std::string(name_suffix)).c_str());
  q2->setRange("q2_range1", 0,    2.95);
  q2->setRange("q2_range2", 2.95, 5.9);
  q2->setRange("q2_range3", 5.9,  8.85);
  q2->setRange("q2_range4", 8.85, 11.8);


  
  //Try to see all the components (why shoudlnt you?)
  
  TString components = "";
  //RooAbsData* channelData1 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  //RooAbsData* channelData2 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  //RooAbsData* channelData3 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  //RooAbsData* channelData4 = (RooAbsData*) data->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel()));
  RooAbsData* channelData1 = (RooAbsData*) data[idx->getLabel()];
  RooAbsData* channelData2 = (RooAbsData*) data[idx->getLabel()];
  RooAbsData* channelData3 = (RooAbsData*) data[idx->getLabel()];
  RooAbsData* channelData4 = (RooAbsData*) data[idx->getLabel()];

  channelData1->plotOn(frame1, RooFit::DrawOption("ZP"), RooFit::CutRange("q2_range1"));
  channelData2->plotOn(frame2, RooFit::DrawOption("ZP"), RooFit::CutRange("q2_range2"));
  channelData3->plotOn(frame3, RooFit::DrawOption("ZP"), RooFit::CutRange("q2_range3"));
  channelData4->plotOn(frame4, RooFit::DrawOption("ZP"), RooFit::CutRange("q2_range4"));
  int colour_count=0;
  
  for (int i = 0; i < active_components_names.size(); ++i)
  {
    std::cout << active_components_names[i] << std::endl;
    if (components != "")
    {
      components+=",";
    }
    components+="*";
    components+=active_components_names[i];
    components+="*";

    int colour_code = component_color(active_components_names[i]);
    
    //model_sumpdf->plotOn(frame1,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData1),RooFit::ProjectionRange("q2_range1"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot1").Data()));     
    //model_sumpdf->plotOn(frame2,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData2),RooFit::ProjectionRange("q2_range2"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot2").Data()));     
    //model_sumpdf->plotOn(frame3,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData3),RooFit::ProjectionRange("q2_range3"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot3").Data()));     
    //model_sumpdf->plotOn(frame4,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjWData(*idx,*channelData4),RooFit::ProjectionRange("q2_range4"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot4").Data()));
    model_sumpdf->plotOn(frame1,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjectionRange("q2_range1"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot1").Data()));     
    model_sumpdf->plotOn(frame2,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjectionRange("q2_range2"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot2").Data()));     
    model_sumpdf->plotOn(frame3,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjectionRange("q2_range3"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot3").Data()));     
    model_sumpdf->plotOn(frame4,RooFit::Slice(*idx, TString(idx->getLabel())),RooFit::ProjectionRange("q2_range4"), RooFit::DrawOption("F"),RooFit::LineColor(colour_code),RooFit::FillColor(colour_code),RooFit::Components(components), RooFit::MoveToBack(), RooFit::Name(active_components_names[i]+TString("_plot4").Data()));    
   
    colour_count++;
  }

  //model_sumpdf->plotOn(frame1, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*channelData1),RooFit::ProjectionRange("q2_range1"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  //model_sumpdf->plotOn(frame2, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*channelData2),RooFit::ProjectionRange("q2_range2"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  //model_sumpdf->plotOn(frame3, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*channelData3),RooFit::ProjectionRange("q2_range3"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  //model_sumpdf->plotOn(frame4, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjWData(*idx,*channelData4),RooFit::ProjectionRange("q2_range4"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  model_sumpdf->plotOn(frame1, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjectionRange("q2_range1"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  model_sumpdf->plotOn(frame2, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjectionRange("q2_range2"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  model_sumpdf->plotOn(frame3, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjectionRange("q2_range3"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));
  model_sumpdf->plotOn(frame4, RooFit::Slice(*idx, TString(idx->getLabel())), RooFit::ProjectionRange("q2_range4"),RooFit::DrawOption("L"),RooFit::LineColor(kBlack), RooFit::LineWidth(1.5));

  //-----Redraw the data points. I don't really understand why I have to do this, to be honest.
  channelData1     = (RooAbsData*) channelData1->reduce(TString(" (obs_x_RDplus_kinematic")+TString(name_suffix)+TString("> 0      & obs_x_RDplus_kinematic")+TString(name_suffix)+TString(" <= 2.95)"));
  channelData2     = (RooAbsData*) channelData2->reduce(TString(" (obs_x_RDplus_kinematic")+TString(name_suffix)+TString("> 2.95   & obs_x_RDplus_kinematic")+TString(name_suffix)+TString(" <= 5.90)"));
  channelData3     = (RooAbsData*) channelData3->reduce(TString(" (obs_x_RDplus_kinematic")+TString(name_suffix)+TString("> 5.90   & obs_x_RDplus_kinematic")+TString(name_suffix)+TString(" <= 8.85)"));
  channelData4     = (RooAbsData*) channelData4->reduce(TString(" (obs_x_RDplus_kinematic")+TString(name_suffix)+TString("> 8.85   & obs_x_RDplus_kinematic")+TString(name_suffix)+TString(" <= 11.80)"));

  channelData1->plotOn(frame1, RooFit::DrawOption("ZP"));
  channelData2->plotOn(frame2, RooFit::DrawOption("ZP"));
  channelData3->plotOn(frame3, RooFit::DrawOption("ZP"));
  channelData4->plotOn(frame4, RooFit::DrawOption("ZP"));
  

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
  
  
  //----Construct the MCstat plot
  std::cout << "-----Creating the MCStatPlot------" << std::endl;   
  MCStatPlot* MCplot1 = new MCStatPlot("MCStatPlot4", model_sumpdf);
  MCStatPlot* MCplot2 = new MCStatPlot("MCStatPlot3", model_sumpdf);
  MCStatPlot* MCplot3 = new MCStatPlot("MCStatPlot2", model_sumpdf);
  MCStatPlot* MCplot4 = new MCStatPlot("MCStatPlot1", model_sumpdf);
  
  MCplot1->SetMap(error_histos);
  MCplot2->SetMap(error_histos);
  MCplot3->SetMap(error_histos);
  MCplot4->SetMap(error_histos);

  MCplot1->SetHistoNameSuffix(TString(name_suffix));
  MCplot2->SetHistoNameSuffix(TString(name_suffix));
  MCplot3->SetHistoNameSuffix(TString(name_suffix));
  MCplot4->SetHistoNameSuffix(TString(name_suffix));
  
  MCplot1->SetNormalizationMode("data");
  MCplot2->SetNormalizationMode("data");
  MCplot3->SetNormalizationMode("data");
  MCplot4->SetNormalizationMode("data");
  
  //TODO: Insert some code in the MCStat plot in order to define some additional cuts to put on the templates and on the data
  MCplot1->SetCut("x", 0.00, 2.95);
  MCplot2->SetCut("x", 2.95, 5.90);
  MCplot3->SetCut("x", 5.90, 8.86);
  MCplot4->SetCut("x", 8.85, 11.8);
  
  TH1F* MCplot_histogram_projected1;
  TH1F* MCplot_histogram_projected2;
  TH1F* MCplot_histogram_projected3;
  TH1F* MCplot_histogram_projected4;

  if(TString(obs->GetName()).Contains("obs_x")){
    MCplot_histogram_projected1 = (TH1F*) MCplot1->GetHistogram1D("x");
    MCplot_histogram_projected2 = (TH1F*) MCplot2->GetHistogram1D("x");
    MCplot_histogram_projected3 = (TH1F*) MCplot3->GetHistogram1D("x");
    MCplot_histogram_projected4 = (TH1F*) MCplot4->GetHistogram1D("x");
  }
  else if(TString(obs->GetName()).Contains("obs_y")){
    MCplot_histogram_projected1 = (TH1F*) MCplot1->GetHistogram1D("y");
    MCplot_histogram_projected2 = (TH1F*) MCplot2->GetHistogram1D("y");
    MCplot_histogram_projected3 = (TH1F*) MCplot3->GetHistogram1D("y");
    MCplot_histogram_projected4 = (TH1F*) MCplot4->GetHistogram1D("y");
  }
  else if (TString(obs->GetName()).Contains("obs_z")){
    MCplot_histogram_projected1 = (TH1F*) MCplot1->GetHistogram1D("z");
    MCplot_histogram_projected2 = (TH1F*) MCplot2->GetHistogram1D("z");
    MCplot_histogram_projected3 = (TH1F*) MCplot3->GetHistogram1D("z");
    MCplot_histogram_projected4 = (TH1F*) MCplot4->GetHistogram1D("z");
  }

  MCplot_histogram_projected1->SetFillStyle(3244);
  MCplot_histogram_projected2->SetFillStyle(3244);
  MCplot_histogram_projected3->SetFillStyle(3244);
  MCplot_histogram_projected4->SetFillStyle(3244);
  
  MCplot_histogram_projected1->SetFillColor(9);
  MCplot_histogram_projected2->SetFillColor(9);
  MCplot_histogram_projected3->SetFillColor(9);
  MCplot_histogram_projected4->SetFillColor(9);
  
  MCplot_histogram_projected1->SetMarkerColor(0);
  MCplot_histogram_projected2->SetMarkerColor(0);
  MCplot_histogram_projected3->SetMarkerColor(0);
  MCplot_histogram_projected4->SetMarkerColor(0);
  
  MCplot_histogram_projected1->SetMarkerSize(0.00001);
  MCplot_histogram_projected2->SetMarkerSize(0.00001);
  MCplot_histogram_projected3->SetMarkerSize(0.00001);
  MCplot_histogram_projected4->SetMarkerSize(0.00001);
  std::cout << "----------------------------------" << std::endl;
  

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

  Double_t xmin = -5.;
  Double_t xmax = 5.;
  TF1* upperLine_function = new TF1("upperLine", "3.", xmin, xmax);
  TF1* mediumLine_function= new TF1("mediumLine", "0.", xmin, xmax);
  TF1* lowerLine_function = new TF1("lowerLine", "-3.",xmin, xmax);
  RooAbsReal* upperLine_roofit = RooFit::bindFunction(upperLine_function, *obs);
  RooAbsReal* mediumLine_roofit= RooFit::bindFunction(mediumLine_function,*obs);
  RooAbsReal* lowerLine_roofit = RooFit::bindFunction(lowerLine_function, *obs);

  upperLine_roofit->plotOn(pframe1, RooFit::LineStyle(3), RooFit::LineColor(12));
  upperLine_roofit->plotOn(pframe2, RooFit::LineStyle(3), RooFit::LineColor(12));
  upperLine_roofit->plotOn(pframe3, RooFit::LineStyle(3), RooFit::LineColor(12));
  upperLine_roofit->plotOn(pframe4, RooFit::LineStyle(3), RooFit::LineColor(12));
  
  mediumLine_roofit->plotOn(pframe1, RooFit::LineStyle(3),RooFit::LineColor(12));
  mediumLine_roofit->plotOn(pframe2, RooFit::LineStyle(3),RooFit::LineColor(12));
  mediumLine_roofit->plotOn(pframe3, RooFit::LineStyle(3),RooFit::LineColor(12));
  mediumLine_roofit->plotOn(pframe4, RooFit::LineStyle(3),RooFit::LineColor(12));
  
  lowerLine_roofit->plotOn(pframe1, RooFit::LineStyle(3), RooFit::LineColor(12));
  lowerLine_roofit->plotOn(pframe2, RooFit::LineStyle(3), RooFit::LineColor(12));
  lowerLine_roofit->plotOn(pframe3, RooFit::LineStyle(3), RooFit::LineColor(12));
  lowerLine_roofit->plotOn(pframe4, RooFit::LineStyle(3), RooFit::LineColor(12));

  c1->cd(2);
  pframe1->Draw();

  c2->cd(2);
  pframe2->Draw();

  c3->cd(2);
  pframe3->Draw();

  c4->cd(2);
  pframe4->Draw();

  //total_templates_uncertainty.Draw("e2,same");
  c1->cd(2);
  MCplot_histogram_projected1->Draw("e2,same");
  c2->cd(2);
  MCplot_histogram_projected2->Draw("e2,same");
  c3->cd(2);
  MCplot_histogram_projected3->Draw("e2,same");
  c4->cd(2);
  MCplot_histogram_projected4->Draw("e2,same");

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
  q2bin1->DrawLatexNDC(0.4, 0.90, "0 < q^{2} < 2.95 GeV^{2}/c^{4}");
  c2->cd(1);
  lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
  q2bin2->DrawLatexNDC(0.4, 0.90, "2.95 < q^{2} < 5.9 GeV^{2}/c^{4}");
  c3->cd(1);
  lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
  q2bin3->DrawLatexNDC(0.4, 0.90, "5.9 < q^{2} < 8.85 GeV^{2}/c^{4}");
  c4->cd(1);
  lhcb->DrawLatexNDC(0.75, 0.90, "LHCb Preliminary");
  q2bin4->DrawLatexNDC(0.4, 0.90, "8.85 < q^{2} < 11.8 GeV^{2}/c^{4}");


  c1->Print(folder+TString("/")+fitvar+TString("_q2Bin1_")+output_file_name+TString(name_suffix)+TString(".pdf"));
  c2->Print(folder+TString("/")+fitvar+TString("_q2Bin2_")+output_file_name+TString(name_suffix)+TString(".pdf"));
  c3->Print(folder+TString("/")+fitvar+TString("_q2Bin3_")+output_file_name+TString(name_suffix)+TString(".pdf"));
  c4->Print(folder+TString("/")+fitvar+TString("_q2Bin4_")+output_file_name+TString(name_suffix)+TString(".pdf"));
  
  c1->Print(folder+TString("/")+fitvar+TString("_q2Bin1_")+output_file_name+TString(name_suffix)+TString(".C"));
  c2->Print(folder+TString("/")+fitvar+TString("_q2Bin2_")+output_file_name+TString(name_suffix)+TString(".C"));
  c3->Print(folder+TString("/")+fitvar+TString("_q2Bin3_")+output_file_name+TString(name_suffix)+TString(".C"));
  c4->Print(folder+TString("/")+fitvar+TString("_q2Bin4_")+output_file_name+TString(name_suffix)+TString(".C"));
  
}

void AddSampleToChannel(RooStats::HistFactory::Channel* channel, RooStats::HistFactory::Sample* sample, std::string isocat_string, std::string sample_name)
{
  /***
  Function used to add the samples to the single channel.à
  It needs in input the name of the sample, the isocat number 
  the sample and the channel pointers.
  It then checks if the sample is one of the banned components in the channel
  and if it is not it adds the sample to the channel
  **/
  
  bool isVetoed;

  isVetoed = IsSampleVetoed(isocat_string, sample_name);

  if (!isVetoed){
    channel->AddSample(*sample);
    std::cout << "Adding sample " << sample->GetName() << " to the Channel " << channel->GetName() << std::endl;
  }
  else{
    std::cout <<"NOT Adding sample " << sample->GetName() << " to the Channel " << channel->GetName() << " since is IS VETOED!" << std::endl; 
  }

  return;
}

int main(int argc, char** argv)
{ 
   
   InputData inputData[] = 
    {
      //{"NPER      ",    "N Events per Run                                             ",35939},
      //{"NRUNS     ",    "N Runs per Job                                               ",500},
      //{"SEED      ",    "random seed                                                  ",0},
      //{"TAG       ",    "tagger,1 means ost, 2 means sst, 3 means combined            ",3},
      //{"NN        ",     "NN cut             ",                                         0.26},
     {"SYST      ",     "systematic error NO             ",                            0},
     {"TOYN      ",     "NO of generated toy             ",                            0},
     //{"SCAN      ",     "2d scan NO             ",                                     0},
     //{"BsPtMin   ",     "Minimal Pt(Bs) for scale factor studies.",                   -99999. },
     //{"BsPtMax   ",     "Maximal Pt(Bs) for scale factor studies.",                   -99999. },
     //{"BsIsoMin  ",     "Minimal Iso(Bs) for scale factor studies.",                  -99999. },
     //{"BsIsoMax  ",     "Maximal Iso(Bs) for scale factor studies.",                  -99999. },
     //{"BsZ0Min  ",      "Minimal Z0 (Bs) for scale factor studies.",                  -99999. },
     //{"BsZ0Max  ",      "Maximal Z0 (Bs) for scale factor studies.",                  -99999. },
     //{"BsDrMin  ",      "Minimal dR (Bs) for scale factor studies.",                  -99999. },
     //{"BsDrMax  ",      "Maximal dR (Bs) for scale factor studies.",                  -99999. },
     //{"SigmaCtCorrLevel", "Level of sigma_ct scale factor correction.",               0 }, 
     //{"NNMax     ",     "Maximal NN cut                           ",                  2.0 },
     //{"NNMin     ",     "Minimal NN cut                           ",                  0.26 },
     //{"POINT     ",     "toy MC number                            ",                  -1 },
     //{"POINT2     ",    "toy MC number                            ",                  -1 },
     //{"LowerSbMin",     "Low  end of the lower sideband.",                             5.28894},
     //{"LowerSbMax",     "High end of the lower sideband.",                             5.31472},
     //{"UpperSbMin",     "Low  end of the upper sideband.",                             5.41784}, 
     //{"UpperSbMax",     "High end of the upper sideband.",                             5.44362}
    };

     //--- Make a new config, then let it parse the command line arguments.
     //
   FitConfig * cfg = FitConfig::instance(inputData, sizeof(inputData)/sizeof(InputData));
   int status = cfg->parseCommandLine( argc, argv );
   if (status != 0) {
     exit (status);    // there was a parsing error, abort!
   }

   int correct_for_sweights = 1;
   std::string year = "2016";
   if (cfg->justplot()) {correct_for_sweights = 0;}
year needs to be changed to be an input option

   std::vector<std::string> plotnames = cfg->names();
   std::vector<std::string> labels    = cfg->labels();

//        po::options_description desc("Allowed options");
//        desc.add_options()
//            ("MCstat","Activate MC stat");
//
//        po::variables_map vm;
//        po::store(po::parse_command_line(ac, av, desc), vm);
//        po::notify(vm);

  //gROOT->ProcessLine(".x ./lhcbStyle.C");
  
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;
  //TFile* q = new TFile("Templates.root"); // This is indeed not needed at this stage!!!
  RooStats::HistFactory::Measurement meas("RDplus","RDplus");
  meas.SetOutputFilePrefix(std::string(cfg->outputFolder()) + std::string("/"));
  meas.SetExportOnly(kTRUE);

  std::vector<std::string> filenames = cfg->parameters(); //Names of the files in which the initial parameters are saved
  std::vector<std::string> isocats   = cfg->isocat(); //List of the categories that are going to be fit
  if (filenames.size() == 0) {
    std::cout << "No parameters input file - taking default" << std::endl;
  }

  //--- Set initial parameters from the file with parameters from each of the fit categories
  std::vector<std::ifstream*> files_vector;
  std::vector<RDplus_parameters*> toyparams_vector;

  std::vector<std::string> name_suffix_vector;
  std::vector<std::string> latex_suffix_vector;
  char name_suffix[256];
  for (unsigned int i=0; i<isocats.size(); ++i){
    if      (strcmp(isocats[i].c_str(),"0")==0) {
      name_suffix_vector.push_back(std::string("_0_FullSample"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{ALL}$"));
    }
    else if (strcmp(isocats[i].c_str(),"1")==0) {
      name_suffix_vector.push_back(std::string("_1_IsolatedSample"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{ISO}$"));
    }
    else if (strcmp(isocats[i].c_str(),"2")==0) {
      name_suffix_vector.push_back(std::string("_2_Dst02Dppim"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{1OS}$"));
    }
    else if (strcmp(isocats[i].c_str(),"3")==0) {
      name_suffix_vector.push_back(std::string("_3_Dstp2Dppi0"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{1NEU}$"));
    }
    else if (strcmp(isocats[i].c_str(),"4")==0) {
      name_suffix_vector.push_back(std::string("_4_Dstp2Dppippim"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{2OS}$"));
    }
    else if (strcmp(isocats[i].c_str(),"5")==0) {
      name_suffix_vector.push_back(std::string("_5_Dstst02Dppimpi0"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{1OS,1NEU}$"));
    }
    else if (strcmp(isocats[i].c_str(),"6")==0) {
      name_suffix_vector.push_back(std::string("_6_Dststp2Dppippimpi0"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{2OS,1NEU}$"));
    }
    else if (strcmp(isocats[i].c_str(),"7")==0) {
      name_suffix_vector.push_back(std::string("_7_DoubleCharm_Kpmum"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{1K}$"));
    }
    else if (strcmp(isocats[i].c_str(),"8")==0) {
      name_suffix_vector.push_back(std::string("_8_DoubleCharm_Kppimmum"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{1K}\\pi$"));
    }
    else if (strcmp(isocats[i].c_str(),"9")==0) {
      name_suffix_vector.push_back(std::string("_9_DoubleCharm_KpKmmum"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{2OSK}$"));
    }
    else if (strcmp(isocats[i].c_str(),"10")==0){
      name_suffix_vector.push_back(std::string("_10_high_B0Mass_Iso"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{high}~B^0$ mass"));
    }
    else if (strcmp(isocats[i].c_str(),"11")==0){
      name_suffix_vector.push_back(std::string("_11_high_B0Mass_Full"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{high}~B^0$ mass, ALL"));
    }
    else if (strcmp(isocats[i].c_str(),"12")==0){
      name_suffix_vector.push_back(std::string("_12_DoubleCharm_inclusive"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{DD}$"));
    }
    else if (strcmp(isocats[i].c_str(), "13")==0){
      name_suffix_vector.push_back(std::string("_13_NormalizationEnriched_Iso"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{NORM~ ISO}$"));
    }
    else if (strcmp(isocats[i].c_str(), "14")==0){
      name_suffix_vector.push_back(std::string("_14_NormalizationEnriched_All"));
      latex_suffix_vector.push_back(std::string("$~ \\mathrm{NORM~ ALL}$"));
    }
  }

  // Reading common parameters
  char common_parameters_file[256];
  strcpy (common_parameters_file ,cfg->common_parameters());
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "             Reading common parameters from file  " << common_parameters_file << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::ifstream* common_parameters_stream = new std::ifstream(common_parameters_file);
  RDplus_common_parameters* common_parameters = new RDplus_common_parameters(*common_parameters_stream);
  
  // Reading non common parameters
  for (unsigned int i=0; i<isocats.size(); ++i){
    files_vector.push_back(new std::ifstream(filenames[i].c_str()));
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "             Reading parameters from file  " << filenames[i].c_str() <<  std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    toyparams_vector.push_back(new RDplus_parameters(*files_vector[i], name_suffix_vector[i].c_str(), latex_suffix_vector[i].c_str()));
  } 
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "            Setting parameters of interest from common ones "                    << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "...Parameters of interest..." << std::endl;
  for (unsigned int k=0; k < common_parameters->nparameters(); k++){
    parameter* param = common_parameters->get_parameter(k);
    if (param->get_POI()){
      std::cout<< "parameter to be made POI "<< param->get_name()<< " start value " << param->get_value() << std::endl;
      meas.SetPOI(param->get_name().c_str());
    }
  }
  std::cout << "....Nuisance parameters...." << std::endl;
  for (unsigned int k=0; k < common_parameters->nparameters(); k++){
    parameter* param = common_parameters->get_parameter(k);
    if (param->get_nuisance()){
      std::cout<< "Nuisance parameter "<< param->get_name()<< " start value " << param->get_value() << std::endl;
    }
  }
  std::cout << "....Parameters that are not going to be fitted...." << std::endl;
  for (unsigned int k=0; k < common_parameters->nparameters(); k++){
    parameter* param = common_parameters->get_parameter(k);
    if (!param->get_is_fitted()){
      std::cout<< "Not Fitted Parameter "<< param->get_name()<< " start value " << param->get_value() << std::endl;
    }
  }
  std::cout << "-------------------------------------------------------------------------------" << std::endl;

  for (unsigned int j=0; j < toyparams_vector.size(); j++){
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "            Setting parameters of interest for category " << isocats[j] << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "...Parameters of interest..." << std::endl;
    for (unsigned int k=0; k < toyparams_vector.at(j)->nparameters(); k++){
      parameter* param = toyparams_vector.at(j)->get_parameter(k);
      if (param->get_POI()){
        std::cout<< "parameter to be made POI "<< param->get_name()<< " start value " << param->get_value() << std::endl;
        meas.SetPOI(param->get_name().c_str());
      }
    }
    std::cout << "....Nuisance parameters...." << std::endl;
    for (unsigned int k=0; k < toyparams_vector.at(j)->nparameters(); k++){
      parameter* param = toyparams_vector.at(j)->get_parameter(k);
      if (param->get_nuisance()){
        std::cout<< "Nuisance parameter "<< param->get_name()<< " start value " << param->get_value() << std::endl;
      }
    }
    std::cout << "....Parameters that are not going to be fitted...." << std::endl;
    for (unsigned int k=0; k < toyparams_vector.at(j)->nparameters(); k++){
      parameter* param = toyparams_vector.at(j)->get_parameter(k);
      if (!param->get_is_fitted()){
        std::cout<< "Not Fitted Parameter "<< param->get_name()<< " start value " << param->get_value() << std::endl;
      }
    }
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
  }

  std::cout << "Parameters defined"<< std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;

  std::cout << std::endl;
  
  meas.SetLumi(1.0);
  


  TFile* fin = new TFile(cfg->templates());

  char RDplus_blinding_string[256] = "StringRDplus";
  strcat(RDplus_blinding_string,cfg->projections()); // This creates a different blinding string for each chosen output name
  char RDst_blinding_string[256] = "StringRDst";
  strcat(RDst_blinding_string,cfg->projections());

  std::string name_blinding_delta_RDplus("blinding_delta_RDplus");
  std::string name_blinding_delta_RDst("blinding_delta_RDst");

  RooRealVar* blinding_delta_RDplus = new RooRealVar(name_blinding_delta_RDplus.c_str(), name_blinding_delta_RDplus.c_str(), 0.);
  RooRealVar* blinding_delta_RDst   = new RooRealVar(name_blinding_delta_RDst.c_str()  , name_blinding_delta_RDst.c_str(),   0.);

  if(cfg->BlindR()) {
    double blinding_scale = 0.3;
    common_parameters->RawRDplus.set_blinding(true, blinding_scale, false, RDplus_blinding_string);
    common_parameters->RawRDst.set_blinding(true, blinding_scale, false, RDst_blinding_string);

    blinding_delta_RDplus->setVal(common_parameters->RawRDplus.get_blinding_delta());
    blinding_delta_RDst->setVal(common_parameters->RawRDst.get_blinding_delta());
  }

  RooRealVar* RawRDplus            = new RooRealVar(common_parameters->RawRDplus.get_name().c_str(), common_parameters->RawRDplus.get_description().c_str(), common_parameters->RawRDplus.get_value()+blinding_delta_RDplus->getVal(),  common_parameters->RawRDplus.get_min()+blinding_delta_RDplus->getVal(), common_parameters->RawRDplus.get_max()+blinding_delta_RDplus->getVal());
  RooRealVar* RawRDst              = new RooRealVar(common_parameters->RawRDst.get_name().c_str(), common_parameters->RawRDst.get_description().c_str(), common_parameters->RawRDst.get_value()+blinding_delta_RDst->getVal(),  common_parameters->RawRDst.get_min()+blinding_delta_RDst->getVal(), common_parameters->RawRDst.get_max()+blinding_delta_RDst->getVal());
  RooFormulaVar* RawRDplus_unblind = new RooFormulaVar((std::string("RawRDplus_unblind")).c_str(),(std::string("RawRDplus_unblind")).c_str(),"@0-@1",RooArgList(*RawRDplus,*blinding_delta_RDplus));
  RooFormulaVar* RawRDst_unblind   = new RooFormulaVar((std::string("RawRDst_unblind")).c_str(),  (std::string("RawRDst_unblind")).c_str(),  "@0-@1",RooArgList(*RawRDst,  *blinding_delta_RDst));



    //---Define the vector of the channels
  std::vector <RooStats::HistFactory::Channel*> channels_vector;

  //--Define the vectors for new parameters to be included in the workspace
  std::vector <RooRealVar*> f_Bu_vector;
  std::vector <RooRealVar*> f_Bu_tauonic_vector;
  std::vector <RooRealVar*> f_DD_Bu_vector;
  std::vector <RooRealVar*> f_DD_Bd_vector;
  std::vector <RooRealVar*> f_DD_tauonic_Bu_vector;
  std::vector <RooRealVar*> f_DD_tauonic_Bd_vector;
  std::vector <RooFormulaVar*> fBd2DD_DDs_formula_vector;
  std::vector <RooFormulaVar*> fBu2DD_DDs_formula_vector;
  std::vector <RooFormulaVar*> fBdDDs_formula_vector;
  std::vector <RooFormulaVar*> fBuDDs_formula_vector;
  std::vector <RooFormulaVar*> fBd2DD_MultiBody_formula_vector;
  std::vector <RooFormulaVar*> fBu2DD_MultiBody_formula_vector;

  std::vector <RooFormulaVar*> f_D0star_Bu_formula_vector;
  std::vector <RooFormulaVar*> f_D0star_Bd_formula_vector;
  std::vector <RooFormulaVar*> f_D2star_Bu_formula_vector;
  std::vector <RooFormulaVar*> f_D2star_Bd_formula_vector;
  std::vector <RooFormulaVar*> f_D1_Bu_formula_vector;
  std::vector <RooFormulaVar*> f_D1_Bd_formula_vector;
  std::vector <RooFormulaVar*> f_D1prime_Bu_formula_vector;
  std::vector <RooFormulaVar*> f_D1prime_Bd_formula_vector;

  //Fractions for the Dstst samples
  RooRealVar* f_broad_Bu         = new RooRealVar((common_parameters->f_broad_Bu.get_name()).c_str(),         (common_parameters->f_broad_Bu.get_description()).c_str(),         common_parameters->f_broad_Bu.get_value(),         common_parameters->f_broad_Bu.get_min(),         common_parameters->f_broad_Bu.get_max());
  RooRealVar* f_broad_Bd         = new RooRealVar((common_parameters->f_broad_Bd.get_name()).c_str(),         (common_parameters->f_broad_Bd.get_description()).c_str(),         common_parameters->f_broad_Bd.get_value(),         common_parameters->f_broad_Bd.get_min(),         common_parameters->f_broad_Bd.get_max());
  RooRealVar* f_D1_narrow_Bu     = new RooRealVar((common_parameters->f_D1_narrow_Bu.get_name()).c_str(),     (common_parameters->f_D1_narrow_Bu.get_description()).c_str(),     common_parameters->f_D1_narrow_Bu.get_value(),     common_parameters->f_D1_narrow_Bu.get_min(),     common_parameters->f_D1_narrow_Bu.get_max());
  RooRealVar* f_D1_narrow_Bd     = new RooRealVar((common_parameters->f_D1_narrow_Bd.get_name()).c_str(),     (common_parameters->f_D1_narrow_Bd.get_description()).c_str(),     common_parameters->f_D1_narrow_Bd.get_value(),     common_parameters->f_D1_narrow_Bd.get_min(),     common_parameters->f_D1_narrow_Bd.get_max());
  RooRealVar* f_D1prime_broad_Bu = new RooRealVar((common_parameters->f_D1prime_broad_Bu.get_name()).c_str(), (common_parameters->f_D1prime_broad_Bu.get_description()).c_str(), common_parameters->f_D1prime_broad_Bu.get_value(), common_parameters->f_D1prime_broad_Bu.get_min(), common_parameters->f_D1prime_broad_Bu.get_max());
  RooRealVar* f_D1prime_broad_Bd = new RooRealVar((common_parameters->f_D1prime_broad_Bd.get_name()).c_str(), (common_parameters->f_D1prime_broad_Bd.get_description()).c_str(), common_parameters->f_D1prime_broad_Bd.get_value(), common_parameters->f_D1prime_broad_Bd.get_min(), common_parameters->f_D1prime_broad_Bd.get_max());

  //--Form-factor parameters for B->D BGL
  RooRealVar *delta_ap0 = new RooRealVar((common_parameters->delta_ap0.get_name()).c_str(), (common_parameters->delta_ap0.get_description()).c_str(), common_parameters->delta_ap0.get_value(), common_parameters->delta_ap0.get_min(), common_parameters->delta_ap0.get_max());
  RooRealVar *delta_ap1 = new RooRealVar((common_parameters->delta_ap1.get_name()).c_str(), (common_parameters->delta_ap1.get_description()).c_str(), common_parameters->delta_ap1.get_value(), common_parameters->delta_ap1.get_min(), common_parameters->delta_ap1.get_max());
  RooRealVar *delta_ap2 = new RooRealVar((common_parameters->delta_ap2.get_name()).c_str(), (common_parameters->delta_ap2.get_description()).c_str(), common_parameters->delta_ap2.get_value(), common_parameters->delta_ap2.get_min(), common_parameters->delta_ap2.get_max());
  RooRealVar *delta_ap3 = new RooRealVar((common_parameters->delta_ap3.get_name()).c_str(), (common_parameters->delta_ap3.get_description()).c_str(), common_parameters->delta_ap3.get_value(), common_parameters->delta_ap3.get_min(), common_parameters->delta_ap3.get_max());
  RooRealVar *delta_a00 = new RooRealVar((common_parameters->delta_a00.get_name()).c_str(), (common_parameters->delta_a00.get_description()).c_str(), common_parameters->delta_a00.get_value(), common_parameters->delta_a00.get_min(), common_parameters->delta_a00.get_max());
  RooRealVar *delta_a01 = new RooRealVar((common_parameters->delta_a01.get_name()).c_str(), (common_parameters->delta_a01.get_description()).c_str(), common_parameters->delta_a01.get_value(), common_parameters->delta_a01.get_min(), common_parameters->delta_a01.get_max());
  RooRealVar *delta_a02 = new RooRealVar((common_parameters->delta_a02.get_name()).c_str(), (common_parameters->delta_a02.get_description()).c_str(), common_parameters->delta_a02.get_value(), common_parameters->delta_a02.get_min(), common_parameters->delta_a02.get_max());
  RooRealVar *delta_a03 = new RooRealVar((common_parameters->delta_a03.get_name()).c_str(), (common_parameters->delta_a03.get_description()).c_str(), common_parameters->delta_a03.get_value(), common_parameters->delta_a03.get_min(), common_parameters->delta_a03.get_max());

  //--Gaussian constraints for B->D BGL parameters (in the diagonal basis)
  RooGaussian *delta_ap0_constraint = new RooGaussian("delta_ap0_constraint","delta_ap0_constraint",*delta_ap0,RooFit::RooConst(0.),RooFit::RooConst(BtoDBGL_FFparams_sigma_array[0]));
  RooGaussian *delta_ap1_constraint = new RooGaussian("delta_ap1_constraint","delta_ap1_constraint",*delta_ap1,RooFit::RooConst(0.),RooFit::RooConst(BtoDBGL_FFparams_sigma_array[1]));
  RooGaussian *delta_ap2_constraint = new RooGaussian("delta_ap2_constraint","delta_ap2_constraint",*delta_ap2,RooFit::RooConst(0.),RooFit::RooConst(BtoDBGL_FFparams_sigma_array[2]));
  RooGaussian *delta_a01_constraint = new RooGaussian("delta_a01_constraint","delta_a01_constraint",*delta_a01,RooFit::RooConst(0.),RooFit::RooConst(BtoDBGL_FFparams_sigma_array[3]));
  RooGaussian *delta_a02_constraint = new RooGaussian("delta_a02_constraint","delta_a02_constraint",*delta_a02,RooFit::RooConst(0.),RooFit::RooConst(BtoDBGL_FFparams_sigma_array[4]));

  //--Form-factor parameters for B->D BGL rotated into the HAMMER base (the usual FF base)
  //--This rotation is only needed when gaussian constraints are applied
  RooFormulaVar *delta_ap0_hammerbase = new RooFormulaVar((common_parameters->delta_ap0.get_name()+"_hammerbase").c_str(), (common_parameters->delta_ap0.get_description()+"_hammerbase").c_str(), formula_delta_ap0, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));
  RooFormulaVar *delta_ap1_hammerbase = new RooFormulaVar((common_parameters->delta_ap1.get_name()+"_hammerbase").c_str(), (common_parameters->delta_ap1.get_description()+"_hammerbase").c_str(), formula_delta_ap1, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));
  RooFormulaVar *delta_ap2_hammerbase = new RooFormulaVar((common_parameters->delta_ap2.get_name()+"_hammerbase").c_str(), (common_parameters->delta_ap2.get_description()+"_hammerbase").c_str(), formula_delta_ap2, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));
  RooFormulaVar *delta_a01_hammerbase = new RooFormulaVar((common_parameters->delta_a01.get_name()+"_hammerbase").c_str(), (common_parameters->delta_a01.get_description()+"_hammerbase").c_str(), formula_delta_a01, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));
  RooFormulaVar *delta_a02_hammerbase = new RooFormulaVar((common_parameters->delta_a02.get_name()+"_hammerbase").c_str(), (common_parameters->delta_a02.get_description()+"_hammerbase").c_str(), formula_delta_a02, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));

  //--Kinematic constraint impossed on delta_a00 if the BtoD BGL parameter constraints are applied
  RooFormulaVar *delta_a00_fixed = new RooFormulaVar((common_parameters->delta_a00.get_name()+"_fixed").c_str(), (common_parameters->delta_a00.get_description()+"_fixed").c_str(), formula_delta_a00_fixed, RooArgList(*delta_ap0,*delta_ap1,*delta_ap2,*delta_a01,*delta_a02));

  //--Form-factor parameters for B->D* BGL
  RooRealVar *delta_a0 = new RooRealVar((common_parameters->delta_a0.get_name()).c_str(), (common_parameters->delta_a0.get_description()).c_str(), common_parameters->delta_a0.get_value(), common_parameters->delta_a0.get_min(), common_parameters->delta_a0.get_max());
  RooRealVar *delta_a1 = new RooRealVar((common_parameters->delta_a1.get_name()).c_str(), (common_parameters->delta_a1.get_description()).c_str(), common_parameters->delta_a1.get_value(), common_parameters->delta_a1.get_min(), common_parameters->delta_a1.get_max());
  RooRealVar *delta_a2 = new RooRealVar((common_parameters->delta_a2.get_name()).c_str(), (common_parameters->delta_a2.get_description()).c_str(), common_parameters->delta_a2.get_value(), common_parameters->delta_a2.get_min(), common_parameters->delta_a2.get_max());
  RooRealVar *delta_b0 = new RooRealVar((common_parameters->delta_b0.get_name()).c_str(), (common_parameters->delta_b0.get_description()).c_str(), common_parameters->delta_b0.get_value(), common_parameters->delta_b0.get_min(), common_parameters->delta_b0.get_max());
  RooRealVar *delta_b1 = new RooRealVar((common_parameters->delta_b1.get_name()).c_str(), (common_parameters->delta_b1.get_description()).c_str(), common_parameters->delta_b1.get_value(), common_parameters->delta_b1.get_min(), common_parameters->delta_b1.get_max());
  RooRealVar *delta_b2 = new RooRealVar((common_parameters->delta_b2.get_name()).c_str(), (common_parameters->delta_b2.get_description()).c_str(), common_parameters->delta_b2.get_value(), common_parameters->delta_b2.get_min(), common_parameters->delta_b2.get_max());
  RooRealVar *delta_c1 = new RooRealVar((common_parameters->delta_c1.get_name()).c_str(), (common_parameters->delta_c1.get_description()).c_str(), common_parameters->delta_c1.get_value(), common_parameters->delta_c1.get_min(), common_parameters->delta_c1.get_max());
  RooRealVar *delta_c2 = new RooRealVar((common_parameters->delta_c2.get_name()).c_str(), (common_parameters->delta_c2.get_description()).c_str(), common_parameters->delta_c2.get_value(), common_parameters->delta_c2.get_min(), common_parameters->delta_c2.get_max());
  RooRealVar *delta_d0 = new RooRealVar((common_parameters->delta_d0.get_name()).c_str(), (common_parameters->delta_d0.get_description()).c_str(), common_parameters->delta_d0.get_value(), common_parameters->delta_d0.get_min(), common_parameters->delta_d0.get_max());
  RooRealVar *delta_d1 = new RooRealVar((common_parameters->delta_d1.get_name()).c_str(), (common_parameters->delta_d1.get_description()).c_str(), common_parameters->delta_d1.get_value(), common_parameters->delta_d1.get_min(), common_parameters->delta_d1.get_max());

  //--Gaussian constraints for B->D* BGL parameters (in the diagonal basis)
  RooGaussian *delta_a0_constraint = new RooGaussian("delta_a0_constraint","delta_a0_constraint",*delta_a0,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[0]));
  RooGaussian *delta_a1_constraint = new RooGaussian("delta_a1_constraint","delta_a1_constraint",*delta_a1,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[1]));
  RooGaussian *delta_b0_constraint = new RooGaussian("delta_b0_constraint","delta_b0_constraint",*delta_b0,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[2]));
  RooGaussian *delta_b1_constraint = new RooGaussian("delta_b1_constraint","delta_b1_constraint",*delta_b1,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[3]));
  RooGaussian *delta_c1_constraint = new RooGaussian("delta_c1_constraint","delta_c1_constraint",*delta_c1,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[4]));
  RooGaussian *delta_c2_constraint = new RooGaussian("delta_c2_constraint","delta_c2_constraint",*delta_c2,RooFit::RooConst(0.),RooFit::RooConst(BtoDstBGL_FFparams_sigma_array[5]));

  //--Form-factor parameters for B->Dst BGL rotated into the HAMMER base (the usual FF base)
  //--This rotation is only needed when gaussian constraints are applied
  RooFormulaVar *delta_a0_hammerbase = new RooFormulaVar((common_parameters->delta_a0.get_name()+"_hammerbase").c_str(), (common_parameters->delta_a0.get_description()+"_hammerbase").c_str(), formula_delta_a0, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));
  RooFormulaVar *delta_a1_hammerbase = new RooFormulaVar((common_parameters->delta_a1.get_name()+"_hammerbase").c_str(), (common_parameters->delta_a1.get_description()+"_hammerbase").c_str(), formula_delta_a1, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));
  RooFormulaVar *delta_b0_hammerbase = new RooFormulaVar((common_parameters->delta_b0.get_name()+"_hammerbase").c_str(), (common_parameters->delta_b0.get_description()+"_hammerbase").c_str(), formula_delta_b0, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));
  RooFormulaVar *delta_b1_hammerbase = new RooFormulaVar((common_parameters->delta_b1.get_name()+"_hammerbase").c_str(), (common_parameters->delta_b1.get_description()+"_hammerbase").c_str(), formula_delta_b1, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));
  RooFormulaVar *delta_c1_hammerbase = new RooFormulaVar((common_parameters->delta_c1.get_name()+"_hammerbase").c_str(), (common_parameters->delta_c1.get_description()+"_hammerbase").c_str(), formula_delta_c1, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));
  RooFormulaVar *delta_c2_hammerbase = new RooFormulaVar((common_parameters->delta_c2.get_name()+"_hammerbase").c_str(), (common_parameters->delta_c2.get_description()+"_hammerbase").c_str(), formula_delta_c2, RooArgList(*delta_a0,*delta_a1,*delta_b0,*delta_b1,*delta_c1,*delta_c2));

  delta_ap0->setConstant(1);
  delta_ap1->setConstant(1);
  delta_ap2->setConstant(1);
  delta_ap3->setConstant(1);
  delta_a00->setConstant(1);
  delta_a01->setConstant(1);
  delta_a02->setConstant(1);
  delta_a03->setConstant(1);
  delta_a0->setConstant(1);
  delta_a1->setConstant(1);
  delta_a2->setConstant(1);
  delta_b0->setConstant(1);
  delta_b1->setConstant(1);
  delta_b2->setConstant(1);
  delta_c1->setConstant(1);
  delta_c2->setConstant(1);
  delta_d0->setConstant(1);
  delta_d1->setConstant(1);

  if (cfg->UseHAMMER()) {
    if (cfg->FloatFFcoefs_BtoD()) {
      if (cfg->ConstrainBtoDBGLParams()) {delta_ap0->setConstant(0);}
      else {delta_ap0->setVal(0.);} // If no external constraint is applied, this parameter is fixed, to account for the reduction of degrees of freedom in the model normalisation (it was found to be the delta parameter with the smallest uncertainty in the fit).
      delta_ap1->setConstant(0);
      if (cfg->ConstrainBtoDBGLParams()) {delta_ap2->setConstant(0);}
      else {delta_ap2->setVal(-BtoDBGL_ap_vector_central_values[2]);} // Reduced number of parameters when no external constraints are applied. Value of the delta set to cancel out the central value set in HAMMER.
      if (!cfg->ConstrainBtoDBGLParams()) {delta_a00->setConstant(0);}
      delta_a01->setConstant(0);
      if (cfg->ConstrainBtoDBGLParams()) {delta_a02->setConstant(0);}
      else {delta_a02->setVal(-BtoDBGL_a0_vector_central_values[2]);} // Reduced number of parameters when no external constraints are applied. Value of the delta set to cancel out the central value set in HAMMER.
    }
    if (cfg->FloatFFcoefs_BtoDst()) {
      delta_a0->setConstant(0);
      delta_a1->setConstant(0);
      delta_b0->setConstant(0);
      delta_b1->setConstant(0);
      delta_c1->setConstant(0);
      delta_c2->setConstant(0);
      delta_d0->setConstant(0);
      delta_d1->setConstant(0);
    }
  }


  //--Configure HAMMER for B02DpMuNu
  std::string WCprocessname_00 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_00 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_00 = RooArgSet();
  RooArgSet _imWCparamlist_00 = RooArgSet();
  std::string FFprocessname_00 = "BtoD";
  std::string FFmodelname_00 = "BGLVar";
  //std::string FFmodelname_errors_00 = "BGL";
  std::vector<std::string>* _FFparamnames_00 = new std::vector<std::string>{"delta_ap0","delta_ap1","delta_ap2","delta_ap3","delta_a00","delta_a01","delta_a02","delta_a03"};
  RooArgSet _FFparamlist_00 = RooArgSet(*delta_ap0,*delta_ap1,*delta_ap2,*delta_ap3,*delta_a00,*delta_a01,*delta_a02,*delta_a03);
  RooArgSet _FFparamlist_00_hammerbase = RooArgSet(*delta_ap0_hammerbase,*delta_ap1_hammerbase,*delta_ap2_hammerbase,*delta_ap3,*delta_a00_fixed,*delta_a01_hammerbase,*delta_a02_hammerbase,*delta_a03);
  std::string histoname_noerrors_00 = "decay_model";
  std::string histoname_witherrors_00 = "decay_model_with_errors";
  std::string schemename_00 = "Scheme1";

  //--Configure HAMMER for B02DpTauNu
  std::string WCprocessname_01 = "BtoCTauNu";
  std::vector<std::string>* _WCparamnames_01 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_01 = RooArgSet();
  RooArgSet _imWCparamlist_01 = RooArgSet();
  std::string FFprocessname_01 = "BtoD";
  std::string FFmodelname_01 = "BGLVar";
  //std::string FFmodelname_errors_01 = "BGL";
  std::vector<std::string>* _FFparamnames_01 = new std::vector<std::string>{"delta_ap0","delta_ap1","delta_ap2","delta_ap3","delta_a00","delta_a01","delta_a02","delta_a03"};
  RooArgSet _FFparamlist_01 = RooArgSet(*delta_ap0,*delta_ap1,*delta_ap2,*delta_ap3,*delta_a00,*delta_a01,*delta_a02,*delta_a03);
  RooArgSet _FFparamlist_01_hammerbase = RooArgSet(*delta_ap0_hammerbase,*delta_ap1_hammerbase,*delta_ap2_hammerbase,*delta_ap3,*delta_a00_fixed,*delta_a01_hammerbase,*delta_a02_hammerbase,*delta_a03);
  std::string histoname_noerrors_01 = "decay_model";
  std::string histoname_witherrors_01 = "decay_model_with_errors";
  std::string schemename_01 = "Scheme1";

  //--Configure HAMMER for B02DstMuNu
  std::string WCprocessname_10 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_10 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_10 = RooArgSet();
  RooArgSet _imWCparamlist_10 = RooArgSet();
  std::string FFprocessname_10 = "BtoD*";
  std::string FFmodelname_10 = "BGLVar";
  //std::string FFmodelname_errors_10 = "BGL";
  std::vector<std::string>* _FFparamnames_10 = new std::vector<std::string>{"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1"};
  RooArgSet _FFparamlist_10 = RooArgSet(*delta_a0,*delta_a1,*delta_a2,*delta_b0,*delta_b1,*delta_b2,*delta_c1,*delta_c2,*delta_d0); // arg list saturated with 9 elements
  _FFparamlist_10.add(*delta_d1);
  RooArgSet _FFparamlist_10_hammerbase = RooArgSet(*delta_a0_hammerbase,*delta_a1_hammerbase,*delta_a2,*delta_b0_hammerbase,*delta_b1_hammerbase,*delta_b2,*delta_c1_hammerbase,*delta_c2_hammerbase,*delta_d0); // arg list saturated with 9 elements
  _FFparamlist_10_hammerbase.add(*delta_d1);
  std::string histoname_noerrors_10 = "decay_model";
  std::string histoname_witherrors_10 = "decay_model_with_errors";
  std::string schemename_10 = "Scheme1";

  //--Configure HAMMER for B02DstTauNu
  std::string WCprocessname_11 = "BtoCTauNu";
  std::vector<std::string>* _WCparamnames_11 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_11 = RooArgSet();
  RooArgSet _imWCparamlist_11 = RooArgSet();
  std::string FFprocessname_11 = "BtoD*";
  std::string FFmodelname_11 = "BGLVar";
  //std::string FFmodelname_errors_11 = "BGL";
  std::vector<std::string>* _FFparamnames_11 = new std::vector<std::string>{"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1"};
  RooArgSet _FFparamlist_11 = RooArgSet(*delta_a0,*delta_a1,*delta_a2,*delta_b0,*delta_b1,*delta_b2,*delta_c1,*delta_c2,*delta_d0);
  _FFparamlist_11.add(*delta_d1);
  RooArgSet _FFparamlist_11_hammerbase = RooArgSet(*delta_a0_hammerbase,*delta_a1_hammerbase,*delta_a2,*delta_b0_hammerbase,*delta_b1_hammerbase,*delta_b2,*delta_c1_hammerbase,*delta_c2_hammerbase,*delta_d0); // arg list saturated with 9 elements
  _FFparamlist_11_hammerbase.add(*delta_d1);
  std::string histoname_noerrors_11 = "decay_model";
  std::string histoname_witherrors_11 = "decay_model_with_errors";
  std::string schemename_11 = "Scheme1";

  //--Configure HAMMER for Bp2Dstst0MuNu
  std::string WCprocessname_20 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_20 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_20 = RooArgSet();
  RooArgSet _imWCparamlist_20 = RooArgSet();
  std::string FFprocessname_20 = "BD**0*";
  std::string FFmodelname_20 = "LLSW";
  //std::string FFmodelname_errors_20 = "LLSW";
  std::vector<std::string>* _FFparamnames_20 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_20 = RooArgSet();
  std::string histoname_noerrors_20 = "decay_model";
  std::string histoname_witherrors_20 = "decay_model_with_errors";
  std::string schemename_20 = "Scheme1";

  //--Configure HAMMER for Bp2Dstst2MuNu
  std::string WCprocessname_30 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_30 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_30 = RooArgSet();
  RooArgSet _imWCparamlist_30 = RooArgSet();
  std::string FFprocessname_30 = "BD**2*";
  std::string FFmodelname_30 = "LLSW";
  //std::string FFmodelname_errors_30 = "LLSW";
  std::vector<std::string>* _FFparamnames_30 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_30 = RooArgSet();
  std::string histoname_noerrors_30 = "decay_model";
  std::string histoname_witherrors_30 = "decay_model_with_errors";
  std::string schemename_30 = "Scheme1";

  //--Configure HAMMER for Bp2Dstst1MuNu
  std::string WCprocessname_40 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_40 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_40 = RooArgSet();
  RooArgSet _imWCparamlist_40 = RooArgSet();
  std::string FFprocessname_40 = "BD**1";
  std::string FFmodelname_40 = "LLSW";
  //std::string FFmodelname_errors_40 = "LLSW";
  std::vector<std::string>* _FFparamnames_40 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_40 = RooArgSet();
  std::string histoname_noerrors_40 = "decay_model";
  std::string histoname_witherrors_40 = "decay_model_with_errors";
  std::string schemename_40 = "Scheme1";

  //--Configure HAMMER for Bp2Dstst1primeMuNu
  std::string WCprocessname_50 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_50 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_50 = RooArgSet();
  RooArgSet _imWCparamlist_50 = RooArgSet();
  std::string FFprocessname_50 = "BD**1*";
  std::string FFmodelname_50 = "LLSW";
  //std::string FFmodelname_errors_50 = "LLSW";
  std::vector<std::string>* _FFparamnames_50 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_50 = RooArgSet();
  std::string histoname_noerrors_50 = "decay_model";
  std::string histoname_witherrors_50 = "decay_model_with_errors";
  std::string schemename_50 = "Scheme1";

  //--Configure HAMMER for B02Dstst0MuNu
  std::string WCprocessname_60 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_60 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_60 = RooArgSet();
  RooArgSet _imWCparamlist_60 = RooArgSet();
  std::string FFprocessname_60 = "BD**0*";
  std::string FFmodelname_60 = "LLSW";
  //std::string FFmodelname_errors_60 = "LLSW";
  std::vector<std::string>* _FFparamnames_60 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_60 = RooArgSet();
  std::string histoname_noerrors_60 = "decay_model";
  std::string histoname_witherrors_60 = "decay_model_with_errors";
  std::string schemename_60 = "Scheme1";

  //--Configure HAMMER for B02Dstst2MuNu
  std::string WCprocessname_70 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_70 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_70 = RooArgSet();
  RooArgSet _imWCparamlist_70 = RooArgSet();
  std::string FFprocessname_70 = "BD**2*";
  std::string FFmodelname_70 = "LLSW";
  //std::string FFmodelname_errors_70 = "LLSW";
  std::vector<std::string>* _FFparamnames_70 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_70 = RooArgSet();
  std::string histoname_noerrors_70 = "decay_model";
  std::string histoname_witherrors_70 = "decay_model_with_errors";
  std::string schemename_70 = "Scheme1";

  //--Configure HAMMER for B02Dstst1MuNu
  std::string WCprocessname_80 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_80 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_80 = RooArgSet();
  RooArgSet _imWCparamlist_80 = RooArgSet();
  std::string FFprocessname_80 = "BD**1";
  std::string FFmodelname_80 = "LLSW";
  //std::string FFmodelname_errors_80 = "LLSW";
  std::vector<std::string>* _FFparamnames_80 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_80 = RooArgSet();
  std::string histoname_noerrors_80 = "decay_model";
  std::string histoname_witherrors_80 = "decay_model_with_errors";
  std::string schemename_80 = "Scheme1";

  //--Configure HAMMER for B02Dstst1primeMuNu
  std::string WCprocessname_90 = "BtoCMuNu";
  std::vector<std::string>* _WCparamnames_90 = new std::vector<std::string>{};
  RooArgSet _reWCparamlist_90 = RooArgSet();
  RooArgSet _imWCparamlist_90 = RooArgSet();
  std::string FFprocessname_90 = "BD**1*";
  std::string FFmodelname_90 = "LLSW";
  //std::string FFmodelname_errors_90 = "LLSW";
  std::vector<std::string>* _FFparamnames_90 = new std::vector<std::string>{};
  RooArgSet _FFparamlist_90 = RooArgSet();
  std::string histoname_noerrors_90 = "decay_model";
  std::string histoname_witherrors_90 = "decay_model_with_errors";
  std::string schemename_90 = "Scheme1";

  //--Specify the HAMMER files and the HAMMER PDFs for each category
  std::vector<std::vector<std::string>> filenames_00;
  std::vector<std::vector<std::string>> filenames_01;
  std::vector<std::vector<std::string>> filenames_10;
  std::vector<std::vector<std::string>> filenames_11;
  std::vector<std::vector<std::string>> filenames_20;
  std::vector<std::vector<std::string>> filenames_30;
  std::vector<std::vector<std::string>> filenames_40;
  std::vector<std::vector<std::string>> filenames_50;
  std::vector<std::vector<std::string>> filenames_60;
  std::vector<std::vector<std::string>> filenames_70;
  std::vector<std::vector<std::string>> filenames_80;
  std::vector<std::vector<std::string>> filenames_90;
  std::vector<std::string> aux_filename_00_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_01_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_10_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_11_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_20_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_30_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_40_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_50_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_60_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_70_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_80_vector = std::vector<std::string> {};
  std::vector<std::string> aux_filename_90_vector = std::vector<std::string> {};

  std::vector<RooHammerModel> cat_HAMMER_PDFs_00 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_01 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_10 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_11 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_20 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_30 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_40 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_50 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_60 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_70 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_80 = std::vector<RooHammerModel> {};
  std::vector<RooHammerModel> cat_HAMMER_PDFs_90 = std::vector<RooHammerModel> {};

  //--Construct one channel for each category

  char datafile_MC[256];
  char datafile_MC_hammer_init[] = "hammerinithistos.root";
  if (correct_for_sweights) {
    TFile* newdatahistosfile = new TFile("NewDataHistograms.root","recreate");
    newdatahistosfile->Close();
    strcpy (datafile_MC,"NewDataHistograms.root");
    //strcpy (datafile_MC_hammer_init,"NewDataHistograms.root");
  }
  else {
    strcpy (datafile_MC,cfg->templates());
    //strcpy (datafile_MC_hammer_init,"hammerinithistos.root");
  }
  std::string histo_component_names[] = {
    "hNDpMuNu",
    "hNDstMuNu",
    "hNBd2DD_DDs",
    "hNBu2DD_DDs",
    "hNBd2DD_MultiBody",
    "hNBu2DD_MultiBody",
    "hNBd2DststmunuHigher",
    "hNMisID",
    "hNWS",
    "hNDpTauNu",
    "hNDstTauNu",
    "hNLbDLc",
    "hNBdDDs",
    "hNBuDDs",
    "hNBd2DD_MultiBody_m1sigma",
    "hNBd2DD_MultiBody_p1sigma",
    "hNBd2DD_MultiBody_m1sigma_quad",
    "hNBd2DD_MultiBody_p1sigma_quad",
    "hNBu2DD_MultiBody_m1sigma",
    "hNBu2DD_MultiBody_p1sigma",
    "hNBu2DD_MultiBody_m1sigma_quad",
    "hNBu2DD_MultiBody_p1sigma_quad",
    "hNBu2Dststmunu_D1prime",
    "hNBu2Dststmunu_D0star",
    "hNBu2Dststmunu_D1",
    "hNBu2Dststmunu_D2star",
    "hNBd2Dststmunu_D1prime",
    "hNBd2Dststmunu_D0star",
    "hNBd2Dststmunu_D1",
    "hNBd2Dststmunu_D2star",
    "hNBd2DststmunuHigher_m1sigma",
    "hNBd2DststmunuHigher_p1sigma"
  };
  //std::string histo_component_names_hammer[] = {
  //  "HAMMER_inithisto_00_cat",
  //  "HAMMER_inithisto_01_cat",
  //  "HAMMER_inithisto_10_cat",
  //  "HAMMER_inithisto_11_cat",
  //  "HAMMER_inithisto_20_cat",
  //  "HAMMER_inithisto_30_cat",
  //  "HAMMER_inithisto_40_cat",
  //  "HAMMER_inithisto_50_cat",
  //  "HAMMER_inithisto_60_cat",
  //  "HAMMER_inithisto_70_cat",
  //  "HAMMER_inithisto_80_cat",
  //  "HAMMER_inithisto_90_cat"
  //};

  std::vector<TH3D*>* sweightcorrhisto = new std::vector<TH3D*>{};
  //std::vector<TH3D*>* sweightcorrhisto_inv = new std::vector<TH3D*>{};
  //std::vector<RooDataHist> sweightcorr_datahist;
  //std::vector<RooHistPdf> sweightcorr_histpdf;
  //std::vector<RooProdPdf> sweightcorr_model_chan;
  //std::vector<RooStats::HistFactory::NormFactor> sweight_normfactor = std::vector<RooStats::HistFactory::NormFactor>{};

  std::vector<std::map<std::string, HistogramErrors>> error_histos;

  for (unsigned int i=0; i<isocats.size(); i++){
    
    std::string channel_name( std::string("RDplus_kinematic")+name_suffix_vector[i]);

    channels_vector.push_back(new RooStats::HistFactory::Channel(channel_name.c_str()));
    channels_vector[i]->SetStatErrorConfig(1e-5, "Poisson");

    error_histos.push_back(std::map<std::string, HistogramErrors>{});

    //-----------------------Data sample------------------------------
    char hData_name[256];
    strcpy (hData_name,"hNData");
    strcat (hData_name,name_suffix_vector[i].c_str());

    //sweight_normfactor.push_back(RooStats::HistFactory::NormFactor());
    //sweight_normfactor.at(i).SetName(("sweight_normfactor_chan"+name_suffix_vector[i]).c_str());

    if (correct_for_sweights) {
        TFile* olddatafile = new TFile(cfg->templates(),"read");
        TFile* newdatafile = new TFile("NewDataHistograms.root","update");
        TH3D* datahisto = (TH3D*) olddatafile->Get(hData_name);
        //double integral_correction = 1./datahisto->Integral();
        sweightcorrhisto->push_back((TH3D*) datahisto->Clone(("th1_sweightcorr"+name_suffix_vector[i]).c_str()));
        sweightcorrhisto->at(i)->SetDirectory(0);
        sweightcorrhisto->at(i)->Reset();
        //sweightcorrhisto_inv->push_back((TH3D*) datahisto->Clone(("th1_sweightcorr_inv"+name_suffix_vector[i]).c_str()));
        //sweightcorrhisto_inv->at(i)->SetDirectory(0);
        //sweightcorrhisto_inv->at(i)->Reset();
        double corr_num, corr_den;
        for (unsigned int xbin_idx=1; xbin_idx<=datahisto->GetXaxis()->GetNbins(); xbin_idx++){
            for (unsigned int ybin_idx=1; ybin_idx<=datahisto->GetYaxis()->GetNbins(); ybin_idx++){
                for (unsigned int zbin_idx=1; zbin_idx<=datahisto->GetZaxis()->GetNbins(); zbin_idx++){
                    corr_num = datahisto->GetBinContent(xbin_idx,ybin_idx,zbin_idx);
                    corr_den = datahisto->GetBinError(xbin_idx,ybin_idx,zbin_idx)*datahisto->GetBinError(xbin_idx,ybin_idx,zbin_idx);
                    if (corr_num <= 0.) {
                        corr_num = 1.;
                        corr_den = 1.;
                    }
                    sweightcorrhisto->at(i)->SetBinContent(xbin_idx,ybin_idx,zbin_idx,corr_num/corr_den);
                    sweightcorrhisto->at(i)->SetBinError(xbin_idx,ybin_idx,zbin_idx,0.);
                    //sweightcorrhisto_inv->at(i)->SetBinContent(xbin_idx,ybin_idx,zbin_idx,corr_den/corr_num);
                }
            }
        }

        datahisto->Sumw2();
        datahisto->Multiply(sweightcorrhisto->at(i));
        //integral_correction *= datahisto->Integral();
        //sweight_normfactor.at(i).SetVal(integral_correction);
        //sweight_normfactor.at(i).SetLow(integral_correction);
        //sweight_normfactor.at(i).SetHigh(integral_correction);
        datahisto->Write();

        //double integral_correction;
        std::vector<TH3D*>* new_histos_1 = new std::vector<TH3D*>{};
        for (unsigned int k=0; k<size(histo_component_names); k++){
            new_histos_1->push_back((TH3D*) olddatafile->Get((histo_component_names[k]+name_suffix_vector[i]).c_str()));
            //integral_correction = 1./new_histos_1->at(k)->Integral();
            new_histos_1->at(k)->Sumw2();
            new_histos_1->at(k)->Multiply(sweightcorrhisto->at(i));
            //integral_correction *= new_histos_1->at(k)->Integral();
            //new_histos_1->at(k)->Scale(integral_correction);
            new_histos_1->at(k)->Write();
        }
        //if (cfg->UseHAMMER()) {
        //    std::vector<TH3D*>* new_histos_2 = new std::vector<TH3D*>{};
        //    TFile* hammerhistosfile = new TFile("hammerinithistos.root","read");
        //    newdatafile->cd();
        //    for (unsigned int k=0; k<size(histo_component_names_hammer); k++){
        //        new_histos_2->push_back((TH3D*) hammerhistosfile->Get((histo_component_names_hammer[k]+isocats[i]).c_str()));
        //        //integral_correction = 1./new_histos_2->at(k)->Integral();
        //        new_histos_2->at(k)->Multiply(sweightcorrhisto->at(i));
        //        //integral_correction *= new_histos_2->at(k)->Integral();
        //        //new_histos_2->at(k)->Scale(integral_correction);
        //        new_histos_2->at(k)->Write();
        //    }
        //    hammerhistosfile->Close();
        //}

        newdatafile->Close();
        olddatafile->Close();

        channels_vector[i]->SetData(hData_name, "NewDataHistograms.root");
    }
    else {
        channels_vector[i]->SetData(hData_name, cfg->templates());
    }

    HistogramErrors Data_Errors;
    Data_Errors.filename = cfg->templates();
    Data_Errors.histogramname = hData_name;

    error_histos[i].emplace(hData_name, Data_Errors);

  }

  if (cfg->UseHAMMER()) {
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "                    Constructing the HAMMER-based PDFs "        << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    for (int icat=0; icat<isocats.size(); icat++){
      aux_filename_00_vector.clear();
      aux_filename_01_vector.clear();
      aux_filename_10_vector.clear();
      aux_filename_11_vector.clear();
      aux_filename_20_vector.clear();
      aux_filename_30_vector.clear();
      aux_filename_40_vector.clear();
      aux_filename_50_vector.clear();
      aux_filename_60_vector.clear();
      aux_filename_70_vector.clear();
      aux_filename_80_vector.clear();
      aux_filename_90_vector.clear();
      for (int ithread=0; ithread<cfg->num_HAMMER_files(); ithread++){
        aux_filename_00_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc0l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_01_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc0l1_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_10_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc1l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_11_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc1l1_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_20_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc2l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_30_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc3l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_40_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc4l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_50_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc5l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_60_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc6l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_70_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc7l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_80_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc8l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
        aux_filename_90_vector.push_back(std::string(cfg->HAMMER_files_dir())+"HammerFile_"+year+"_Hc9l0_cat"+isocats[icat]+"_thread"+std::to_string(ithread)+".root");
      }
      filenames_00.push_back(std::vector<std::string>(aux_filename_00_vector));
      filenames_01.push_back(std::vector<std::string>(aux_filename_01_vector));
      filenames_10.push_back(std::vector<std::string>(aux_filename_10_vector));
      filenames_11.push_back(std::vector<std::string>(aux_filename_11_vector));
      filenames_20.push_back(std::vector<std::string>(aux_filename_20_vector));
      filenames_30.push_back(std::vector<std::string>(aux_filename_30_vector));
      filenames_40.push_back(std::vector<std::string>(aux_filename_40_vector));
      filenames_50.push_back(std::vector<std::string>(aux_filename_50_vector));
      filenames_60.push_back(std::vector<std::string>(aux_filename_60_vector));
      filenames_70.push_back(std::vector<std::string>(aux_filename_70_vector));
      filenames_80.push_back(std::vector<std::string>(aux_filename_80_vector));
      filenames_90.push_back(std::vector<std::string>(aux_filename_90_vector));

      if (cfg->FloatFFcoefs_BtoD() and cfg->ConstrainBtoDBGLParams()) {
        cat_HAMMER_PDFs_00.push_back(RooHammerModel(((std::string) "HAMMER_PDF_00_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_00_cat"+isocats[icat]).c_str(),WCprocessname_00,_WCparamnames_00,_reWCparamlist_00,_imWCparamlist_00,FFprocessname_00,FFmodelname_00,_FFparamnames_00,_FFparamlist_00_hammerbase,&filenames_00.at(icat),histoname_noerrors_00,histoname_witherrors_00,schemename_00));
        cat_HAMMER_PDFs_01.push_back(RooHammerModel(((std::string) "HAMMER_PDF_01_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_01_cat"+isocats[icat]).c_str(),WCprocessname_01,_WCparamnames_01,_reWCparamlist_01,_imWCparamlist_01,FFprocessname_01,FFmodelname_01,_FFparamnames_01,_FFparamlist_01_hammerbase,&filenames_01.at(icat),histoname_noerrors_01,histoname_witherrors_01,schemename_01));
      }
      else {
        cat_HAMMER_PDFs_00.push_back(RooHammerModel(((std::string) "HAMMER_PDF_00_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_00_cat"+isocats[icat]).c_str(),WCprocessname_00,_WCparamnames_00,_reWCparamlist_00,_imWCparamlist_00,FFprocessname_00,FFmodelname_00,_FFparamnames_00,_FFparamlist_00,&filenames_00.at(icat),histoname_noerrors_00,histoname_witherrors_00,schemename_00));
        cat_HAMMER_PDFs_01.push_back(RooHammerModel(((std::string) "HAMMER_PDF_01_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_01_cat"+isocats[icat]).c_str(),WCprocessname_01,_WCparamnames_01,_reWCparamlist_01,_imWCparamlist_01,FFprocessname_01,FFmodelname_01,_FFparamnames_01,_FFparamlist_01,&filenames_01.at(icat),histoname_noerrors_01,histoname_witherrors_01,schemename_01));
      }
      if (cfg->FloatFFcoefs_BtoDst() and cfg->ConstrainBtoDstBGLParams()) {
        cat_HAMMER_PDFs_10.push_back(RooHammerModel(((std::string) "HAMMER_PDF_10_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_10_cat"+isocats[icat]).c_str(),WCprocessname_10,_WCparamnames_10,_reWCparamlist_10,_imWCparamlist_10,FFprocessname_10,FFmodelname_10,_FFparamnames_10,_FFparamlist_10_hammerbase,&filenames_10.at(icat),histoname_noerrors_10,histoname_witherrors_10,schemename_10));
        cat_HAMMER_PDFs_11.push_back(RooHammerModel(((std::string) "HAMMER_PDF_11_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_11_cat"+isocats[icat]).c_str(),WCprocessname_11,_WCparamnames_11,_reWCparamlist_11,_imWCparamlist_11,FFprocessname_11,FFmodelname_11,_FFparamnames_11,_FFparamlist_11_hammerbase,&filenames_11.at(icat),histoname_noerrors_11,histoname_witherrors_11,schemename_11));
      }
      else {
        cat_HAMMER_PDFs_10.push_back(RooHammerModel(((std::string) "HAMMER_PDF_10_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_10_cat"+isocats[icat]).c_str(),WCprocessname_10,_WCparamnames_10,_reWCparamlist_10,_imWCparamlist_10,FFprocessname_10,FFmodelname_10,_FFparamnames_10,_FFparamlist_10,&filenames_10.at(icat),histoname_noerrors_10,histoname_witherrors_10,schemename_10));
        cat_HAMMER_PDFs_11.push_back(RooHammerModel(((std::string) "HAMMER_PDF_11_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_11_cat"+isocats[icat]).c_str(),WCprocessname_11,_WCparamnames_11,_reWCparamlist_11,_imWCparamlist_11,FFprocessname_11,FFmodelname_11,_FFparamnames_11,_FFparamlist_11,&filenames_11.at(icat),histoname_noerrors_11,histoname_witherrors_11,schemename_11));
      }
      cat_HAMMER_PDFs_20.push_back(RooHammerModel(((std::string) "HAMMER_PDF_20_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_20_cat"+isocats[icat]).c_str(),WCprocessname_20,_WCparamnames_20,_reWCparamlist_20,_imWCparamlist_20,FFprocessname_20,FFmodelname_20,_FFparamnames_20,_FFparamlist_20,&filenames_20.at(icat),histoname_noerrors_20,histoname_witherrors_20,schemename_20));
      cat_HAMMER_PDFs_30.push_back(RooHammerModel(((std::string) "HAMMER_PDF_30_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_30_cat"+isocats[icat]).c_str(),WCprocessname_30,_WCparamnames_30,_reWCparamlist_30,_imWCparamlist_30,FFprocessname_30,FFmodelname_30,_FFparamnames_30,_FFparamlist_30,&filenames_30.at(icat),histoname_noerrors_30,histoname_witherrors_30,schemename_30));
      cat_HAMMER_PDFs_40.push_back(RooHammerModel(((std::string) "HAMMER_PDF_40_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_40_cat"+isocats[icat]).c_str(),WCprocessname_40,_WCparamnames_40,_reWCparamlist_40,_imWCparamlist_40,FFprocessname_40,FFmodelname_40,_FFparamnames_40,_FFparamlist_40,&filenames_40.at(icat),histoname_noerrors_40,histoname_witherrors_40,schemename_40));
      cat_HAMMER_PDFs_50.push_back(RooHammerModel(((std::string) "HAMMER_PDF_50_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_50_cat"+isocats[icat]).c_str(),WCprocessname_50,_WCparamnames_50,_reWCparamlist_50,_imWCparamlist_50,FFprocessname_50,FFmodelname_50,_FFparamnames_50,_FFparamlist_50,&filenames_50.at(icat),histoname_noerrors_50,histoname_witherrors_50,schemename_50));
      cat_HAMMER_PDFs_60.push_back(RooHammerModel(((std::string) "HAMMER_PDF_60_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_60_cat"+isocats[icat]).c_str(),WCprocessname_60,_WCparamnames_60,_reWCparamlist_60,_imWCparamlist_60,FFprocessname_60,FFmodelname_60,_FFparamnames_60,_FFparamlist_60,&filenames_60.at(icat),histoname_noerrors_60,histoname_witherrors_60,schemename_60));
      cat_HAMMER_PDFs_70.push_back(RooHammerModel(((std::string) "HAMMER_PDF_70_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_70_cat"+isocats[icat]).c_str(),WCprocessname_70,_WCparamnames_70,_reWCparamlist_70,_imWCparamlist_70,FFprocessname_70,FFmodelname_70,_FFparamnames_70,_FFparamlist_70,&filenames_70.at(icat),histoname_noerrors_70,histoname_witherrors_70,schemename_70));
      cat_HAMMER_PDFs_80.push_back(RooHammerModel(((std::string) "HAMMER_PDF_80_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_80_cat"+isocats[icat]).c_str(),WCprocessname_80,_WCparamnames_80,_reWCparamlist_80,_imWCparamlist_80,FFprocessname_80,FFmodelname_80,_FFparamnames_80,_FFparamlist_80,&filenames_80.at(icat),histoname_noerrors_80,histoname_witherrors_80,schemename_80));
      cat_HAMMER_PDFs_90.push_back(RooHammerModel(((std::string) "HAMMER_PDF_90_cat"+isocats[icat]).c_str(),((std::string) "HAMMER_PDF_90_cat"+isocats[icat]).c_str(),WCprocessname_90,_WCparamnames_90,_reWCparamlist_90,_imWCparamlist_90,FFprocessname_90,FFmodelname_90,_FFparamnames_90,_FFparamlist_90,&filenames_90.at(icat),histoname_noerrors_90,histoname_witherrors_90,schemename_90)); 
    }
    std::cout << "-------------------------------------------------------------------------------" << std::endl;  
    std::cout << std::endl;
  }

  //--Correct for sWeights if required
  if (cfg->UseHAMMER() and correct_for_sweights) {
   for (int icat=0; icat<isocats.size(); icat++){
      cat_HAMMER_PDFs_00[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_01[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_10[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_11[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_20[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_30[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_40[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_50[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_60[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_70[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_80[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
      cat_HAMMER_PDFs_90[icat].SetExternalMultiplicativeHistogram(sweightcorrhisto->at(icat));
    }
  }

  //--Store the inital histograms provided by HAMMER in a file
  TFile* fhammerinithistos;
  std::vector<TH3D>* vhammerinithistos;
  if (cfg->UseHAMMER()) {
    fhammerinithistos = new TFile("hammerinithistos.root","recreate");
    vhammerinithistos = new std::vector<TH3D>{};
    for (int icat=0; icat<isocats.size(); icat++){
      std::cout << "FFParamList after constructing " << _FFparamlist_00 << std::endl;
      std::cout << "                               " << _FFparamlist_20 << std::endl;
      vhammerinithistos->push_back(cat_HAMMER_PDFs_00[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_00_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_01[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_01_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_10[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_10_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_11[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_11_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_20[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_20_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_30[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_30_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_40[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_40_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_50[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_50_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_60[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_60_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_70[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_70_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_80[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_80_cat"+isocats[icat]).c_str()));
      vhammerinithistos->push_back(cat_HAMMER_PDFs_90[icat].getHistogramWithErrors(((std::string) "HAMMER_inithisto_90_cat"+isocats[icat]).c_str()));
    }
    for (int ihisto=0; ihisto<vhammerinithistos->size(); ihisto++){vhammerinithistos->at(ihisto).Write();}
    fhammerinithistos->Close();
  }

  for (unsigned int i=0; i<isocats.size(); i++){

    //-----------------------B02DpMuNu------------------------------
    char hNDpMuNu_name[256];
    strcpy (hNDpMuNu_name,"hNDpMuNu");
    strcat (hNDpMuNu_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* DplusMuNu;
    if (cfg->UseHAMMER()) {
        DplusMuNu = new RooStats::HistFactory::Sample(hNDpMuNu_name,"HAMMER_inithisto_00_cat"+isocats[i],datafile_MC_hammer_init);
        
        HistogramErrors DplusMuNu_Errors;
        DplusMuNu_Errors.filename      = "hammerinithistos.root";
        DplusMuNu_Errors.histogramname = "HAMMER_inithisto_00_cat"+isocats[i];
        error_histos[i].emplace(hNDpMuNu_name, DplusMuNu_Errors);
      }
    else {
      TH1D* hDpMu = (TH1D*)fin->Get(hNDpMuNu_name);
      DplusMuNu = new RooStats::HistFactory::Sample(hNDpMuNu_name,hNDpMuNu_name,datafile_MC);
      DplusMuNu->AddNormFactor((std::string("MCDplusMu")+name_suffix_vector[i]).c_str(),1.0/hDpMu->Integral(),0.0,1.0);

      HistogramErrors DplusMuNu_Errors;
      DplusMuNu_Errors.filename      = cfg->templates();
      DplusMuNu_Errors.histogramname = hNDpMuNu_name;
      error_histos[i].emplace(hNDpMuNu_name, DplusMuNu_Errors);

    }
    DplusMuNu->SetNormalizeByTheory(kFALSE);
    DplusMuNu->AddNormFactor(common_parameters->NDpMuNu.get_name(), common_parameters->NDpMuNu.get_value(),  common_parameters->NDpMuNu.get_min(), common_parameters->NDpMuNu.get_max());
    DplusMuNu->AddNormFactor(toyparams_vector[i]->TF_NDpmunu.get_name(), toyparams_vector[i]->TF_NDpmunu.get_value(), toyparams_vector[i]->TF_NDpmunu.get_min(), toyparams_vector[i]->TF_NDpmunu.get_max());
    //if (correct_for_sweights) {DplusMuNu->AddNormFactor(sweight_normfactor.at(i));}
    
    //-----------------------B02DstMuNu------------------------------
    char hNDstMuNu_name[256];
    strcpy (hNDstMuNu_name,"hNDstMuNu");
    strcat (hNDstMuNu_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* DstMuNu;
    if (cfg->UseHAMMER()) {
      DstMuNu = new RooStats::HistFactory::Sample(hNDstMuNu_name,"HAMMER_inithisto_10_cat"+isocats[i],datafile_MC_hammer_init);

      HistogramErrors DstMuNu_Errors;
      DstMuNu_Errors.filename = "hammerinithistos.root";
      DstMuNu_Errors.histogramname = "HAMMER_inithisto_10_cat"+isocats[i];

      error_histos[i].emplace(hNDstMuNu_name, DstMuNu_Errors);
    }
    else {
      TH1D* hDstMu = (TH1D*)fin->Get(hNDstMuNu_name);
      DstMuNu = new RooStats::HistFactory::Sample(hNDstMuNu_name,hNDstMuNu_name,datafile_MC);
      DstMuNu->AddNormFactor((std::string("MCDstMu")+name_suffix_vector[i]).c_str(),1.0/hDstMu->Integral(),0.0,1.0);

      HistogramErrors DstMuNu_Errors;
      DstMuNu_Errors.filename = cfg->templates();
      DstMuNu_Errors.histogramname = hNDstMuNu_name;
      error_histos[i].emplace(hNDstMuNu_name, DstMuNu_Errors); 
    }
    DstMuNu->SetNormalizeByTheory(kFALSE);
    DstMuNu->AddNormFactor(common_parameters->NDstMuNu.get_name(), common_parameters->NDstMuNu.get_value(),  common_parameters->NDstMuNu.get_min(), common_parameters->NDstMuNu.get_max());
    DstMuNu->AddNormFactor(toyparams_vector[i]->TF_NDstmunu.get_name(), toyparams_vector[i]->TF_NDstmunu.get_value(), toyparams_vector[i]->TF_NDstmunu.get_min(), toyparams_vector[i]->TF_NDstmunu.get_max());
    //if (correct_for_sweights) {DstMuNu->AddNormFactor(sweight_normfactor.at(i));}

    //---------------------B02DpTauNu-------------------------------
    char hNDpTauNu_name[256];
    strcpy (hNDpTauNu_name,"hNDpTauNu");
    strcat (hNDpTauNu_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* DplusTauNu;
    if (cfg->UseHAMMER()) {
      DplusTauNu = new RooStats::HistFactory::Sample(hNDpTauNu_name,"HAMMER_inithisto_01_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors DplusTauNu_Errors;
      DplusTauNu_Errors.filename = "hammerinithistos.root";
      DplusTauNu_Errors.histogramname = "HAMMER_inithisto_01_cat"+isocats[i];

      error_histos[i].emplace(hNDpTauNu_name, DplusTauNu_Errors);
    }
    
    else {
      TH1D* hDpTau = (TH1D*)fin->Get(hNDpTauNu_name);
      DplusTauNu = new RooStats::HistFactory::Sample(hNDpTauNu_name,hNDpTauNu_name,datafile_MC);
      DplusTauNu->AddNormFactor((std::string("MCDplusTau")+name_suffix_vector[i]).c_str(),1.0/hDpTau->Integral(),0.0,1.0);

      HistogramErrors DplusTauNu_Errors;
      DplusTauNu_Errors.filename = cfg->templates();
      DplusTauNu_Errors.histogramname = hNDpTauNu_name;
      error_histos[i].emplace(hNDpTauNu_name, DplusTauNu_Errors); 
    }
    DplusTauNu->SetNormalizeByTheory(kFALSE);
    DplusTauNu->AddNormFactor(common_parameters->NDpMuNu.get_name(), common_parameters->NDpMuNu.get_value(),  common_parameters->NDpMuNu.get_min(), common_parameters->NDpMuNu.get_max());
    DplusTauNu->AddNormFactor(toyparams_vector[i]->TF_NDptaunu.get_name(), toyparams_vector[i]->TF_NDptaunu.get_value(), toyparams_vector[i]->TF_NDptaunu.get_min(), toyparams_vector[i]->TF_NDptaunu.get_max());
    DplusTauNu->AddNormFactor((std::string("RawRDplus_dummy")).c_str(), 0.04, 0., 1.);
    //if (correct_for_sweights) {DplusTauNu->AddNormFactor(sweight_normfactor.at(i));}

    //--------------------B02DstTauNu--------------------------------
    char hNDstTauNu_name[256];
    strcpy (hNDstTauNu_name,"hNDstTauNu");
    strcat (hNDstTauNu_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* DstTauNu;
    if (cfg->UseHAMMER()) {
      DstTauNu = new RooStats::HistFactory::Sample(hNDstTauNu_name,"HAMMER_inithisto_11_cat"+isocats[i],datafile_MC_hammer_init);
      
      HistogramErrors DstTauNu_Errors;
      DstTauNu_Errors.filename = "hammerinithistos.root";
      DstTauNu_Errors.histogramname = "HAMMER_inithisto_11_cat"+isocats[i];

      error_histos[i].emplace(hNDstTauNu_name, DstTauNu_Errors);
    }
    
    else {
      TH1D* hDstTau = (TH1D*)fin->Get(hNDstTauNu_name);
      DstTauNu = new RooStats::HistFactory::Sample(hNDstTauNu_name,hNDstTauNu_name,datafile_MC);
      DstTauNu->AddNormFactor((std::string("MCDstTau")+name_suffix_vector[i]).c_str(),1.0/hDstTau->Integral(),0.0,1.0);

      HistogramErrors DstTauNu_Errors;
      DstTauNu_Errors.filename = cfg->templates();
      DstTauNu_Errors.histogramname = hNDstTauNu_name;

      error_histos[i].emplace(hNDstTauNu_name, DstTauNu_Errors);

    } 
    DstTauNu->SetNormalizeByTheory(kFALSE);
    DstTauNu->AddNormFactor(common_parameters->NDstMuNu.get_name(), common_parameters->NDstMuNu.get_value(),  common_parameters->NDstMuNu.get_min(), common_parameters->NDstMuNu.get_max());
    DstTauNu->AddNormFactor(toyparams_vector[i]->TF_NDsttaunu.get_name(), toyparams_vector[i]->TF_NDsttaunu.get_value(), toyparams_vector[i]->TF_NDsttaunu.get_min(), toyparams_vector[i]->TF_NDsttaunu.get_max());
    DstTauNu->AddNormFactor((std::string("RawRDst_dummy")).c_str(), 0.04, 0., 1.);
    //if (correct_for_sweights) {DstTauNu->AddNormFactor(sweight_normfactor.at(i));}

    //--------------------DD samples

    RooStats::HistFactory::NormFactor* NDD = new RooStats::HistFactory::NormFactor(); //Common normalization to all the DD samples
    NDD->SetName(toyparams_vector[i]->NDD.get_name());
    NDD->SetVal(toyparams_vector[i]->NDD.get_value());
    NDD->SetLow(toyparams_vector[i]->NDD.get_min());
    NDD->SetHigh(toyparams_vector[i]->NDD.get_max());

    RooStats::HistFactory::NormFactor* NDD_tau = new RooStats::HistFactory::NormFactor();
    NDD_tau->SetName(toyparams_vector[i]->NDD_tau.get_name());
    NDD_tau->SetVal(toyparams_vector[i]->NDD_tau.get_value());
    NDD_tau->SetLow(toyparams_vector[i]->NDD_tau.get_min());
    NDD_tau->SetHigh(toyparams_vector[i]->NDD_tau.get_max());

    RooRealVar* f_Bu            = new RooRealVar(toyparams_vector[i]->f_Bu.get_name().c_str(),    toyparams_vector[i]->f_Bu.get_description().c_str(),    toyparams_vector[i]->f_Bu.get_value(),    toyparams_vector[i]->f_Bu.get_min(),    toyparams_vector[i]->f_Bu.get_max());
    RooRealVar* f_Bu_tauonic    = new RooRealVar((std::string("f_Bu_tauonic")+name_suffix_vector[i]).c_str(),    (std::string("f_Bu_tauonic")+name_suffix_vector[i]).c_str(),    toyparams_vector[i]->f_Bu.get_value(),    toyparams_vector[i]->f_Bu.get_min(),    toyparams_vector[i]->f_Bu.get_max());
    RooRealVar* f_DD_Bu = new RooRealVar(toyparams_vector[i]->f_DD_Bu.get_name().c_str(), toyparams_vector[i]->f_DD_Bu.get_description().c_str(), toyparams_vector[i]->f_DD_Bu.get_value(), toyparams_vector[i]->f_DD_Bu.get_min(), toyparams_vector[i]->f_DD_Bu.get_max()); 
    RooRealVar* f_DD_Bd = new RooRealVar(toyparams_vector[i]->f_DD_Bd.get_name().c_str(), toyparams_vector[i]->f_DD_Bd.get_description().c_str(), toyparams_vector[i]->f_DD_Bd.get_value(), toyparams_vector[i]->f_DD_Bd.get_min(), toyparams_vector[i]->f_DD_Bd.get_max());
    RooRealVar* f_DD_tauonic_Bu = new RooRealVar(toyparams_vector[i]->f_DD_tauonic_Bu.get_name().c_str(), toyparams_vector[i]->f_DD_tauonic_Bu.get_description().c_str(), toyparams_vector[i]->f_DD_tauonic_Bu.get_value(), toyparams_vector[i]->f_DD_tauonic_Bu.get_min(), toyparams_vector[i]->f_DD_tauonic_Bu.get_max()); 
    RooRealVar* f_DD_tauonic_Bd = new RooRealVar(toyparams_vector[i]->f_DD_tauonic_Bd.get_name().c_str(), toyparams_vector[i]->f_DD_tauonic_Bd.get_description().c_str(), toyparams_vector[i]->f_DD_tauonic_Bd.get_value(), toyparams_vector[i]->f_DD_tauonic_Bd.get_min(), toyparams_vector[i]->f_DD_tauonic_Bd.get_max());
    f_Bu_vector.push_back(f_Bu);
    f_Bu_tauonic_vector.push_back(f_Bu_tauonic);
    f_DD_Bu_vector.push_back(f_DD_Bu);
    f_DD_Bd_vector.push_back(f_DD_Bd);
    f_DD_tauonic_Bu_vector.push_back(f_DD_tauonic_Bu);
    f_DD_tauonic_Bd_vector.push_back(f_DD_tauonic_Bd);
    
    //-------------------B02DD_TwoBody_muonic--------------------------------
    char hNBd2DD_DDs_name[256];
    strcpy (hNBd2DD_DDs_name,"hNBd2DD_DDs");
    strcat (hNBd2DD_DDs_name,name_suffix_vector[i].c_str());
    TH1D* hBd2DD_DDs = (TH1D*)fin->Get(hNBd2DD_DDs_name);
    RooStats::HistFactory::Sample* Bd2DD_DDs = new RooStats::HistFactory::Sample(hNBd2DD_DDs_name,hNBd2DD_DDs_name,datafile_MC);
    Bd2DD_DDs->SetNormalizeByTheory(kFALSE);
    Bd2DD_DDs->AddNormFactor((std::string("MCBd2DD_DDs")+name_suffix_vector[i]).c_str(),1.0/hBd2DD_DDs->Integral(),0.0,1.0);
    //Bd2DD_DDs.AddNormFactor(toyparams.NBd2DD_DDs.get_name(),toyparams.NBd2DD_DDs.get_value(),toyparams.NBd2DD_DDs.get_min(),toyparams.NBd2DD_DDs.get_max());
    Bd2DD_DDs->AddNormFactor(*NDD);
    Bd2DD_DDs->AddNormFactor((std::string("fBd2DD_DDs")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {Bd2DD_DDs->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* fBd2DD_DDs_formula = new RooFormulaVar((std::string("fBd2DD_DDs_formula")+name_suffix_vector[i]).c_str(),(std::string("fBd2DD_DDs_formula")+name_suffix_vector[i]).c_str(),"(1.-@0)*@1",RooArgList(*f_Bu_vector[i],*f_DD_Bd_vector[i]));
    fBd2DD_DDs_formula_vector.push_back(fBd2DD_DDs_formula);

    HistogramErrors Bd2DD_DDs_Errors;
    Bd2DD_DDs_Errors.filename = cfg->templates();
    Bd2DD_DDs_Errors.histogramname = hNBd2DD_DDs_name;
    error_histos[i].emplace(hNBd2DD_DDs_name, Bd2DD_DDs_Errors); 
    
    //-------------------Bplus2DD_TwoBody_muonic------------------------------
    char hNBu2DD_DDs_name[256];
    strcpy (hNBu2DD_DDs_name,"hNBu2DD_DDs");
    strcat (hNBu2DD_DDs_name,name_suffix_vector[i].c_str());
    TH1D* hBu2DD_DDs = (TH1D*)fin->Get(hNBu2DD_DDs_name);
    RooStats::HistFactory::Sample* Bu2DD_DDs = new RooStats::HistFactory::Sample(hNBu2DD_DDs_name,hNBu2DD_DDs_name,datafile_MC);
    Bu2DD_DDs->SetNormalizeByTheory(kFALSE);
    Bu2DD_DDs->AddNormFactor((std::string("MCBu2DD_DDs")+name_suffix_vector[i]).c_str(),1.0/hBu2DD_DDs->Integral(),0.0,1.0);
    //Bu2DD_DDs.AddNormFactor(toyparams.NBu2DD_DDs.get_name(),toyparams.NBu2DD_DDs.get_value(),toyparams.NBu2DD_DDs.get_min(),toyparams.NBu2DD_DDs.get_max());
    Bu2DD_DDs->AddNormFactor(*NDD);
    Bu2DD_DDs->AddNormFactor((std::string("fBu2DD_DDs")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {Bu2DD_DDs->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* fBu2DD_DDs_formula = new RooFormulaVar((std::string("fBu2DD_DDs_formula")+name_suffix_vector[i]).c_str(),(std::string("fBu2DD_DDs_formula")+name_suffix_vector[i]).c_str(),"@0*@1",RooArgList(*f_Bu_vector[i],*f_DD_Bu_vector[i]));
    fBu2DD_DDs_formula_vector.push_back(fBu2DD_DDs_formula);

    HistogramErrors Bu2DD_DDs_Errors;
    Bu2DD_DDs_Errors.filename = cfg->templates();
    Bu2DD_DDs_Errors.histogramname = hNBu2DD_DDs_name;
    error_histos[i].emplace(hNBu2DD_DDs_name, Bu2DD_DDs_Errors); 
    
    //-------------------B02DD_tauonic--------------------------------
    char hNBdDDs_name[256];
    strcpy (hNBdDDs_name,"hNBdDDs");
    strcat (hNBdDDs_name,name_suffix_vector[i].c_str());
    TH1D* hBdDDs = (TH1D*)fin->Get(hNBdDDs_name);
    RooStats::HistFactory::Sample* BdDDs = new RooStats::HistFactory::Sample(hNBdDDs_name,hNBdDDs_name,datafile_MC);
    BdDDs->SetNormalizeByTheory(kFALSE);
    BdDDs->AddNormFactor((std::string("MCBdDDs")+name_suffix_vector[i]).c_str(),1.0/hBdDDs->Integral(),0.0,1.0);
    BdDDs->AddNormFactor(*NDD_tau);
    BdDDs->AddNormFactor((std::string("fBdDDs")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {BdDDs->AddNormFactor(sweight_normfactor.at(i));}
    //RooFormulaVar* fBdDDs_formula = new RooFormulaVar((std::string("fBdDDs_formula")+name_suffix_vector[i]).c_str(),(std::string("fBdDDs_formula")+name_suffix_vector[i]).c_str(),"(1.-@0)*@1*@2",RooArgList(*f_Bu_vector[i],*f_DD_Bd_vector[i],*f_DD_tauonic_Bd_vector[i]));
    RooFormulaVar* fBdDDs_formula = new RooFormulaVar((std::string("fBdDDs_formula")+name_suffix_vector[i]).c_str(),(std::string("fBdDDs_formula")+name_suffix_vector[i]).c_str(),"(1.-@0)",RooArgList(*f_Bu_vector[i]));
    fBdDDs_formula_vector.push_back(fBdDDs_formula);

    HistogramErrors BdDDs_Errors;
    BdDDs_Errors.filename = cfg->templates();
    BdDDs_Errors.histogramname = hNBdDDs_name;
    error_histos[i].emplace(hNBdDDs_name, BdDDs_Errors); 
    
    //-------------------Bplus2DD_tauonic------------------------------
    char hNBuDDs_name[256];
    strcpy (hNBuDDs_name,"hNBuDDs");
    strcat (hNBuDDs_name,name_suffix_vector[i].c_str());
    TH1D* hBuDDs = (TH1D*)fin->Get(hNBuDDs_name);
    RooStats::HistFactory::Sample* BuDDs = new RooStats::HistFactory::Sample(hNBuDDs_name,hNBuDDs_name,datafile_MC);
    BuDDs->SetNormalizeByTheory(kFALSE);
    BuDDs->AddNormFactor((std::string("MCBuDDs")+name_suffix_vector[i]).c_str(),1.0/hBuDDs->Integral(),0.0,1.0);
    BuDDs->AddNormFactor(*NDD_tau);
    BuDDs->AddNormFactor((std::string("fBuDDs")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {BuDDs->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* fBuDDs_formula = new RooFormulaVar((std::string("fBuDDs_formula")+name_suffix_vector[i]).c_str(),(std::string("fBuDDs_formula")+name_suffix_vector[i]).c_str(),"@0",RooArgList(*f_Bu_vector[i]));
    fBuDDs_formula_vector.push_back(fBuDDs_formula);

    HistogramErrors in_Errors;
    in_Errors.filename = cfg->templates();
    in_Errors.histogramname = hNBuDDs_name;
    error_histos[i].emplace(hNBuDDs_name, in_Errors); 

    //------------------B02DD_MultiBody----------------------------------
    char hNBd2DD_MultiBody_name[256];
    strcpy (hNBd2DD_MultiBody_name,"hNBd2DD_MultiBody");
    strcat (hNBd2DD_MultiBody_name,name_suffix_vector[i].c_str());
    TH1D* hBd2DD_MultiBody = (TH1D*)fin->Get(hNBd2DD_MultiBody_name);
    RooStats::HistFactory::Sample* Bd2DD_MultiBody = new RooStats::HistFactory::Sample(hNBd2DD_MultiBody_name,hNBd2DD_MultiBody_name,datafile_MC);
    Bd2DD_MultiBody->SetNormalizeByTheory(kFALSE);
    Bd2DD_MultiBody->AddNormFactor((std::string("MCBd2DD_MultiBody")+name_suffix_vector[i]).c_str(),1.0/hBd2DD_MultiBody->Integral(),0.0,1.0);
    //Bd2DD_MultiBody.AddNormFactor(toyparams.NBd2DD_MultiBody.get_name(),toyparams.NBd2DD_MultiBody.get_value(),toyparams.NBd2DD_MultiBody.get_min(),toyparams.NBd2DD_MultiBody.get_max());
    Bd2DD_MultiBody->AddNormFactor(*NDD);
    Bd2DD_MultiBody->AddNormFactor((std::string("fBd2DD_MultiBody")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {Bd2DD_MultiBody->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* fBd2DD_MultiBody_formula = new RooFormulaVar((std::string("fBd2DD_MultiBody_formula")+name_suffix_vector[i]).c_str(),(std::string("fBd2DD_MultiBody_formula")+name_suffix_vector[i]).c_str(),"(1.-@0)*(1.-@1)",RooArgList(*f_Bu_vector[i],*f_DD_Bd_vector[i]));
    fBd2DD_MultiBody_formula_vector.push_back(fBd2DD_MultiBody_formula);

    HistogramErrors Bd2DD_MultiBody_Errors;
    Bd2DD_MultiBody_Errors.filename = cfg->templates();
    Bd2DD_MultiBody_Errors.histogramname = hNBd2DD_MultiBody_name;

    error_histos[i].emplace(hNBd2DD_MultiBody_name, Bd2DD_MultiBody_Errors); 

    if(cfg->varyshapes())
    {
      char hNBd2DD_MultiBody_m1sigma_name[256];
      strcpy (hNBd2DD_MultiBody_m1sigma_name, "hNBd2DD_MultiBody_m1sigma");
      strcat (hNBd2DD_MultiBody_m1sigma_name,name_suffix_vector[i].c_str());
      char hNBd2DD_MultiBody_p1sigma_name[256];
      strcpy (hNBd2DD_MultiBody_p1sigma_name, "hNBd2DD_MultiBody_p1sigma");
      strcat (hNBd2DD_MultiBody_p1sigma_name,name_suffix_vector[i].c_str());
      //Bd2DD_MultiBody->AddHistoSys((std::string("hNBd2DD_MultiBody_variation")+name_suffix_vector[i]).c_str(),hNBd2DD_MultiBody_m1sigma_name,cfg->templates(),"",hNBd2DD_MultiBody_p1sigma_name,cfg->templates(),"");
      Bd2DD_MultiBody->AddHistoSys("MultiBody_variation",hNBd2DD_MultiBody_m1sigma_name,cfg->templates(),"",hNBd2DD_MultiBody_p1sigma_name,cfg->templates(),"");
      
      char hNBd2DD_MultiBody_m1sigma_quad_name[256];
      strcpy (hNBd2DD_MultiBody_m1sigma_quad_name, "hNBd2DD_MultiBody_m1sigma_quad");
      strcat (hNBd2DD_MultiBody_m1sigma_quad_name,name_suffix_vector[i].c_str());
      char hNBd2DD_MultiBody_p1sigma_quad_name[256];
      strcpy (hNBd2DD_MultiBody_p1sigma_quad_name, "hNBd2DD_MultiBody_p1sigma_quad");
      strcat (hNBd2DD_MultiBody_p1sigma_quad_name,name_suffix_vector[i].c_str());
      //Bd2DD_MultiBody->AddHistoSys((std::string("hNBd2DD_MultiBody_quad_variation")+name_suffix_vector[i]).c_str(),hNBd2DD_MultiBody_m1sigma_quad_name,cfg->templates(),"",hNBd2DD_MultiBody_p1sigma_quad_name,cfg->templates(),"");
      Bd2DD_MultiBody->AddHistoSys("MultiBody_quad_variation",hNBd2DD_MultiBody_m1sigma_quad_name,cfg->templates(),"",hNBd2DD_MultiBody_p1sigma_quad_name,cfg->templates(),"");
    }    
    
    //------------------Bu2DD_MultiBody----------------------------------
    char hNBu2DD_MultiBody_name[256];
    strcpy (hNBu2DD_MultiBody_name,"hNBu2DD_MultiBody");
    strcat (hNBu2DD_MultiBody_name,name_suffix_vector[i].c_str());
    TH1D* hBu2DD_MultiBody = (TH1D*)fin->Get(hNBu2DD_MultiBody_name);
    RooStats::HistFactory::Sample* Bu2DD_MultiBody = new RooStats::HistFactory::Sample(hNBu2DD_MultiBody_name,hNBu2DD_MultiBody_name,datafile_MC);
    Bu2DD_MultiBody->SetNormalizeByTheory(kFALSE);
    Bu2DD_MultiBody->AddNormFactor((std::string("MCBu2DD_MultiBody")+name_suffix_vector[i]).c_str(),1.0/hBu2DD_MultiBody->Integral(),0.0,1.0);
    //Bu2DD_MultiBody.AddNormFactor(toyparams.NBu2DD_MultiBody.get_name(),toyparams.NBu2DD_MultiBody.get_value(),toyparams.NBu2DD_MultiBody.get_min(),toyparams.NBu2DD_MultiBody.get_max());
    Bu2DD_MultiBody->AddNormFactor(*NDD);
    Bu2DD_MultiBody->AddNormFactor((std::string("fBu2DD_MultiBody")+name_suffix_vector[i]).c_str(),0.5,0.0,1.0);
    //if (correct_for_sweights) {Bu2DD_MultiBody->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* fBu2DD_MultiBody_formula = new RooFormulaVar((std::string("fBu2DD_MultiBody_formula")+name_suffix_vector[i]).c_str(),(std::string("fBu2DD_MultiBody_formula")+name_suffix_vector[i]).c_str(),"@0*(1.-@1)",RooArgList(*f_Bu_vector[i],*f_DD_Bu_vector[i]));
    fBu2DD_MultiBody_formula_vector.push_back(fBu2DD_MultiBody_formula);

    HistogramErrors Bu2DD_MultiBody_Errors;
    Bu2DD_MultiBody_Errors.filename = cfg->templates();
    Bu2DD_MultiBody_Errors.histogramname = hNBu2DD_MultiBody_name;

    error_histos[i].emplace(hNBu2DD_MultiBody_name, Bu2DD_MultiBody_Errors); 

    if(cfg->varyshapes())
    {
      char hNBu2DD_MultiBody_m1sigma_name[256];
      strcpy (hNBu2DD_MultiBody_m1sigma_name, "hNBu2DD_MultiBody_m1sigma");
      strcat (hNBu2DD_MultiBody_m1sigma_name,name_suffix_vector[i].c_str());
      char hNBu2DD_MultiBody_p1sigma_name[256];
      strcpy (hNBu2DD_MultiBody_p1sigma_name, "hNBu2DD_MultiBody_p1sigma");
      strcat (hNBu2DD_MultiBody_p1sigma_name,name_suffix_vector[i].c_str());
      //Bu2DD_MultiBody->AddHistoSys((std::string("hNBu2DD_MultiBody_variation")+name_suffix_vector[i]).c_str(),hNBu2DD_MultiBody_m1sigma_name,cfg->templates(),"",hNBu2DD_MultiBody_p1sigma_name,cfg->templates(),"");
      Bu2DD_MultiBody->AddHistoSys("MultiBody_variation",hNBu2DD_MultiBody_m1sigma_name,cfg->templates(),"",hNBu2DD_MultiBody_p1sigma_name,cfg->templates(),"");
      
      char hNBu2DD_MultiBody_m1sigma_quad_name[256];
      strcpy (hNBu2DD_MultiBody_m1sigma_quad_name, "hNBu2DD_MultiBody_m1sigma_quad");
      strcat (hNBu2DD_MultiBody_m1sigma_quad_name,name_suffix_vector[i].c_str());
      char hNBu2DD_MultiBody_p1sigma_quad_name[256];
      strcpy (hNBu2DD_MultiBody_p1sigma_quad_name, "hNBu2DD_MultiBody_p1sigma_quad");
      strcat (hNBu2DD_MultiBody_p1sigma_quad_name,name_suffix_vector[i].c_str());
      //Bu2DD_MultiBody->AddHistoSys((std::string("hNBu2DD_MultiBody_quad_variation")+name_suffix_vector[i]).c_str(),hNBu2DD_MultiBody_m1sigma_quad_name,cfg->templates(),"",hNBu2DD_MultiBody_p1sigma_quad_name,cfg->templates(),"");
      Bu2DD_MultiBody->AddHistoSys("MultiBody_quad_variation",hNBu2DD_MultiBody_m1sigma_quad_name,cfg->templates(),"",hNBu2DD_MultiBody_p1sigma_quad_name,cfg->templates(),"");
      //TODO: Check if giving the same name to two different HistSys variations is okay
    }
    //--------------------------Bplus2Dstst---------------------------
    char hNBu2Dststmunu_D1prime_name[256];
    strcpy (hNBu2Dststmunu_D1prime_name,"hNBu2Dststmunu_D1prime");
    strcat (hNBu2Dststmunu_D1prime_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bu2Dstst_D1prime;
    if (cfg->UseHAMMER()) {
      Bu2Dstst_D1prime = new RooStats::HistFactory::Sample(hNBu2Dststmunu_D1prime_name,"HAMMER_inithisto_50_cat"+isocats[i],datafile_MC_hammer_init);

      HistogramErrors Bu2Dstst_D1prime_Errors;
      Bu2Dstst_D1prime_Errors.filename = "hammerinithistos.root";
      Bu2Dstst_D1prime_Errors.histogramname = "HAMMER_inithisto_50_cat"+isocats[i];

      error_histos[i].emplace(hNBu2Dststmunu_D1prime_name, Bu2Dstst_D1prime_Errors);
    }
    else {
      TH1D* hBu2Dstst_D1prime = (TH1D*)fin->Get(hNBu2Dststmunu_D1prime_name);
      Bu2Dstst_D1prime = new RooStats::HistFactory::Sample (hNBu2Dststmunu_D1prime_name,hNBu2Dststmunu_D1prime_name,datafile_MC);
      Bu2Dstst_D1prime->AddNormFactor((std::string("MCBu2Dstst_D1prime")+name_suffix_vector[i]).c_str(),1.0/hBu2Dstst_D1prime->Integral(),0.0,1.0);

      HistogramErrors Bu2Dstst_D1prime_Errors;
      Bu2Dstst_D1prime_Errors.filename = cfg->templates();
      Bu2Dstst_D1prime_Errors.histogramname = hNBu2Dststmunu_D1prime_name;

      error_histos[i].emplace(hNBu2Dststmunu_D1prime_name, Bu2Dstst_D1prime_Errors);
    }
    Bu2Dstst_D1prime->SetNormalizeByTheory(kFALSE);
    Bu2Dstst_D1prime->AddNormFactor(toyparams_vector[i]->NBu2Dstst.get_name(),toyparams_vector[i]->NBu2Dstst.get_value(),toyparams_vector[i]->NBu2Dstst.get_min(),toyparams_vector[i]->NBu2Dstst.get_max());
    Bu2Dstst_D1prime->AddNormFactor((std::string("f_broad_D1prime_Bu")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bu2Dstst_D1prime->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bu_D1prime.get_name(), toyparams_vector[i]->TF_Dstst_Bu_D1prime.get_value(), toyparams_vector[i]->TF_Dstst_Bu_D1prime.get_min(), toyparams_vector[i]->TF_Dstst_Bu_D1prime.get_max());
    //if (correct_for_sweights) {Bu2Dstst_D1prime->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D1prime_Bu_formula = new RooFormulaVar((std::string("f_D1prime_Bu_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D1prime_Bu_formula")+name_suffix_vector[i]).c_str(), "@0*@1", RooArgList(*f_broad_Bu, *f_D1prime_broad_Bu));
    f_D1prime_Bu_formula_vector.push_back(f_D1prime_Bu_formula);

    char hNBu2Dststmunu_D0star_name[256];
    strcpy (hNBu2Dststmunu_D0star_name,"hNBu2Dststmunu_D0star");
    strcat (hNBu2Dststmunu_D0star_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bu2Dstst_D0star;
    if (cfg->UseHAMMER()) {
      Bu2Dstst_D0star = new RooStats::HistFactory::Sample(hNBu2Dststmunu_D0star_name,"HAMMER_inithisto_20_cat"+isocats[i],datafile_MC_hammer_init);

      HistogramErrors Bu2Dstst_D0star_Errors;
      Bu2Dstst_D0star_Errors.filename = "hammerinithistos.root";
      Bu2Dstst_D0star_Errors.histogramname = "HAMMER_inithisto_20_cat"+isocats[i];

      error_histos[i].emplace(hNBu2Dststmunu_D0star_name, Bu2Dstst_D0star_Errors);
    }
    else {
      TH1D* hBu2Dstst_D0star = (TH1D*)fin->Get(hNBu2Dststmunu_D0star_name);
      Bu2Dstst_D0star = new RooStats::HistFactory::Sample (hNBu2Dststmunu_D0star_name,hNBu2Dststmunu_D0star_name,datafile_MC);
      Bu2Dstst_D0star->AddNormFactor((std::string("MCBu2Dstst_D0star")+name_suffix_vector[i]).c_str(),1.0/hBu2Dstst_D0star->Integral(),0.0,1.0);

      HistogramErrors Bu2Dstst_D0star_Errors;
      Bu2Dstst_D0star_Errors.filename = cfg->templates();
      Bu2Dstst_D0star_Errors.histogramname = hNBu2Dststmunu_D0star_name;

      error_histos[i].emplace(hNBu2Dststmunu_D0star_name, Bu2Dstst_D0star_Errors);
    }
    Bu2Dstst_D0star->SetNormalizeByTheory(kFALSE);
    Bu2Dstst_D0star->AddNormFactor(toyparams_vector[i]->NBu2Dstst.get_name(),toyparams_vector[i]->NBu2Dstst.get_value(),toyparams_vector[i]->NBu2Dstst.get_min(),toyparams_vector[i]->NBu2Dstst.get_max());
    Bu2Dstst_D0star->AddNormFactor((std::string("f_broad_D0star_Bu")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bu2Dstst_D0star->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bu_D0.get_name(), toyparams_vector[i]->TF_Dstst_Bu_D0.get_value(), toyparams_vector[i]->TF_Dstst_Bu_D0.get_min(), toyparams_vector[i]->TF_Dstst_Bu_D0.get_max());
    //if (correct_for_sweights) {Bu2Dstst_D0star->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D0star_Bu_formula = new RooFormulaVar((std::string("f_D0star_Bu_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D0star_Bu_formula")+name_suffix_vector[i]).c_str(), "@0*(1.-@1)", RooArgList(*f_broad_Bu, *f_D1prime_broad_Bu));  
    f_D0star_Bu_formula_vector.push_back(f_D0star_Bu_formula);

    char hNBu2Dststmunu_D1_name[256];
    strcpy (hNBu2Dststmunu_D1_name,"hNBu2Dststmunu_D1");
    strcat (hNBu2Dststmunu_D1_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bu2Dstst_D1;
    if (cfg->UseHAMMER()) {
      Bu2Dstst_D1 = new RooStats::HistFactory::Sample(hNBu2Dststmunu_D1_name,"HAMMER_inithisto_40_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bu2Dstst_D1_Errors;
      Bu2Dstst_D1_Errors.filename = "hammerinithistos.root";
      Bu2Dstst_D1_Errors.histogramname = "HAMMER_inithisto_40_cat"+isocats[i];

      error_histos[i].emplace(hNBu2Dststmunu_D1_name, Bu2Dstst_D1_Errors);
    }
    else {
      TH1D* hBu2Dstst_D1 = (TH1D*)fin->Get(hNBu2Dststmunu_D1_name);
      Bu2Dstst_D1 = new RooStats::HistFactory::Sample (hNBu2Dststmunu_D1_name,hNBu2Dststmunu_D1_name,datafile_MC);
      Bu2Dstst_D1->AddNormFactor((std::string("MCBu2Dstst_D1")+name_suffix_vector[i]).c_str(),1.0/hBu2Dstst_D1->Integral(),0.0,1.0);

      HistogramErrors Bu2Dstst_D1_Errors;
      Bu2Dstst_D1_Errors.filename = cfg->templates();
      Bu2Dstst_D1_Errors.histogramname = hNBu2Dststmunu_D1_name;

      error_histos[i].emplace(hNBu2Dststmunu_D1_name, Bu2Dstst_D1_Errors);
    }
    Bu2Dstst_D1->SetNormalizeByTheory(kFALSE);
    Bu2Dstst_D1->AddNormFactor(toyparams_vector[i]->NBu2Dstst.get_name(),toyparams_vector[i]->NBu2Dstst.get_value(),toyparams_vector[i]->NBu2Dstst.get_min(),toyparams_vector[i]->NBu2Dstst.get_max());
    Bu2Dstst_D1->AddNormFactor((std::string("f_narrow_D1_Bu")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bu2Dstst_D1->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bu_D1.get_name(),toyparams_vector[i]->TF_Dstst_Bu_D1.get_value(), toyparams_vector[i]->TF_Dstst_Bu_D1.get_min(), toyparams_vector[i]->TF_Dstst_Bu_D1.get_max());
    //if (correct_for_sweights) {Bu2Dstst_D1->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D1_Bu_formula = new RooFormulaVar((std::string("f_D1_Bu_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D1_Bu_formula")+name_suffix_vector[i]).c_str(), "(1.-@0)*@1", RooArgList(*f_broad_Bu, *f_D1_narrow_Bu));
    f_D1_Bu_formula_vector.push_back(f_D1_Bu_formula);

    char hNBu2Dststmunu_D2star_name[256];
    strcpy (hNBu2Dststmunu_D2star_name,"hNBu2Dststmunu_D2star");
    strcat (hNBu2Dststmunu_D2star_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bu2Dstst_D2star;
    if (cfg->UseHAMMER()) {
      Bu2Dstst_D2star = new RooStats::HistFactory::Sample(hNBu2Dststmunu_D2star_name,"HAMMER_inithisto_30_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bu2Dstst_D2star_Errors;
      Bu2Dstst_D2star_Errors.filename = "hammerinithistos.root";
      Bu2Dstst_D2star_Errors.histogramname = "HAMMER_inithisto_30_cat"+isocats[i];

      error_histos[i].emplace(hNBu2Dststmunu_D2star_name,Bu2Dstst_D2star_Errors); 
    }
    else {
      TH1D* hBu2Dstst_D2star = (TH1D*)fin->Get(hNBu2Dststmunu_D2star_name);
      Bu2Dstst_D2star = new RooStats::HistFactory::Sample (hNBu2Dststmunu_D2star_name,hNBu2Dststmunu_D2star_name,datafile_MC);
      Bu2Dstst_D2star->AddNormFactor((std::string("MCBu2Dstst_D2star")+name_suffix_vector[i]).c_str(),1.0/hBu2Dstst_D2star->Integral(),0.0,1.0);
    
      HistogramErrors Bu2Dstst_D2star_Errors;
      Bu2Dstst_D2star_Errors.filename = cfg->templates();
      Bu2Dstst_D2star_Errors.histogramname = hNBu2Dststmunu_D2star_name;

      error_histos[i].emplace(hNBu2Dststmunu_D2star_name,Bu2Dstst_D2star_Errors);
    }
    Bu2Dstst_D2star->SetNormalizeByTheory(kFALSE);
    Bu2Dstst_D2star->AddNormFactor(toyparams_vector[i]->NBu2Dstst.get_name(),toyparams_vector[i]->NBu2Dstst.get_value(),toyparams_vector[i]->NBu2Dstst.get_min(),toyparams_vector[i]->NBu2Dstst.get_max());
    Bu2Dstst_D2star->AddNormFactor((std::string("f_narrow_D2star_Bu")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bu2Dstst_D2star->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bu_D2.get_name(),toyparams_vector[i]->TF_Dstst_Bu_D2.get_value(),toyparams_vector[i]->TF_Dstst_Bu_D2.get_min(), toyparams_vector[i]->TF_Dstst_Bu_D2.get_max());
    //if (correct_for_sweights) {Bu2Dstst_D2star->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D2star_Bu_formula = new RooFormulaVar((std::string("f_D2star_Bu_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D2star_Bu_formula")+name_suffix_vector[i]).c_str(), "(1.-@0)*(1.-@1)", RooArgList(*f_broad_Bu, *f_D1_narrow_Bu));
    f_D2star_Bu_formula_vector.push_back(f_D2star_Bu_formula);

    //-------------------------B02Dstst-------------------------------

    char hNBd2Dststmunu_D1prime_name[256];
    strcpy (hNBd2Dststmunu_D1prime_name,"hNBd2Dststmunu_D1prime");
    strcat (hNBd2Dststmunu_D1prime_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bd2Dstst_D1prime;
    if (cfg->UseHAMMER()) {
      Bd2Dstst_D1prime = new RooStats::HistFactory::Sample(hNBd2Dststmunu_D1prime_name,"HAMMER_inithisto_90_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bd2Dstst_D1prime_Errors;
      Bd2Dstst_D1prime_Errors.filename = "hammerinithistos.root";
      Bd2Dstst_D1prime_Errors.histogramname = "HAMMER_inithisto_90_cat"+isocats[i];

      error_histos[i].emplace(hNBd2Dststmunu_D1prime_name,Bd2Dstst_D1prime_Errors); 
    }
    else {
      TH1D* hBd2Dstst_D1prime = (TH1D*)fin->Get(hNBd2Dststmunu_D1prime_name);
      Bd2Dstst_D1prime = new RooStats::HistFactory::Sample (hNBd2Dststmunu_D1prime_name,hNBd2Dststmunu_D1prime_name,datafile_MC);
      Bd2Dstst_D1prime->AddNormFactor((std::string("MCBd2Dstst_D1prime")+name_suffix_vector[i]).c_str(),1.0/hBd2Dstst_D1prime->Integral(),0.0,1.0);
    
      HistogramErrors Bd2Dstst_D1prime_Errors;
      Bd2Dstst_D1prime_Errors.filename = cfg->templates();
      Bd2Dstst_D1prime_Errors.histogramname = hNBd2Dststmunu_D1prime_name;

      error_histos[i].emplace(hNBd2Dststmunu_D1prime_name,Bd2Dstst_D1prime_Errors);
    }
    Bd2Dstst_D1prime->SetNormalizeByTheory(kFALSE);
    Bd2Dstst_D1prime->AddNormFactor(toyparams_vector[i]->NBd2Dstst.get_name(),toyparams_vector[i]->NBd2Dstst.get_value(),toyparams_vector[i]->NBd2Dstst.get_min(),toyparams_vector[i]->NBd2Dstst.get_max());
    Bd2Dstst_D1prime->AddNormFactor((std::string("f_broad_D1prime_Bd")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bd2Dstst_D1prime->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bd_D1prime.get_name(), toyparams_vector[i]->TF_Dstst_Bd_D1prime.get_value(), toyparams_vector[i]->TF_Dstst_Bd_D1prime.get_min(), toyparams_vector[i]->TF_Dstst_Bd_D1prime.get_max());
    //if (correct_for_sweights) {Bd2Dstst_D1prime->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D1prime_Bd_formula = new RooFormulaVar((std::string("f_D1prime_Bd_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D1prime_Bd_formula")+name_suffix_vector[i]).c_str(), "@0*@1", RooArgList(*f_broad_Bd, *f_D1prime_broad_Bd));
    f_D1prime_Bd_formula_vector.push_back(f_D1prime_Bd_formula);

    char hNBd2Dststmunu_D0star_name[256];
    strcpy (hNBd2Dststmunu_D0star_name,"hNBd2Dststmunu_D0star");
    strcat (hNBd2Dststmunu_D0star_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bd2Dstst_D0star;
    if (cfg->UseHAMMER()) {
      Bd2Dstst_D0star = new RooStats::HistFactory::Sample(hNBd2Dststmunu_D0star_name,"HAMMER_inithisto_60_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bd2Dstst_D0star_Errors;
      Bd2Dstst_D0star_Errors.filename = "hammerinithistos.root";
      Bd2Dstst_D0star_Errors.histogramname = "HAMMER_inithisto_60_cat"+isocats[i];

      error_histos[i].emplace(hNBd2Dststmunu_D0star_name,Bd2Dstst_D0star_Errors); 
    }
    else {
      TH1D* hBd2Dstst_D0star = (TH1D*)fin->Get(hNBd2Dststmunu_D0star_name);
      Bd2Dstst_D0star = new RooStats::HistFactory::Sample (hNBd2Dststmunu_D0star_name,hNBd2Dststmunu_D0star_name,datafile_MC);
      Bd2Dstst_D0star->AddNormFactor((std::string("MCBd2Dstst_D0star")+name_suffix_vector[i]).c_str(),1.0/hBd2Dstst_D0star->Integral(),0.0,1.0);
    
      HistogramErrors Bd2Dstst_D0star_Errors;
      Bd2Dstst_D0star_Errors.filename = cfg->templates();
      Bd2Dstst_D0star_Errors.histogramname = hNBd2Dststmunu_D0star_name;

      error_histos[i].emplace(hNBd2Dststmunu_D0star_name,Bd2Dstst_D0star_Errors);
    }
    Bd2Dstst_D0star->SetNormalizeByTheory(kFALSE);
    Bd2Dstst_D0star->AddNormFactor(toyparams_vector[i]->NBd2Dstst.get_name(),toyparams_vector[i]->NBd2Dstst.get_value(),toyparams_vector[i]->NBd2Dstst.get_min(),toyparams_vector[i]->NBd2Dstst.get_max());
    Bd2Dstst_D0star->AddNormFactor((std::string("f_broad_D0star_Bd")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bd2Dstst_D0star->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bd_D0.get_name(), toyparams_vector[i]->TF_Dstst_Bd_D0.get_value(), toyparams_vector[i]->TF_Dstst_Bd_D0.get_min(), toyparams_vector[i]->TF_Dstst_Bd_D0.get_max());
    //if (correct_for_sweights) {Bd2Dstst_D0star->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D0star_Bd_formula = new RooFormulaVar((std::string("f_D0star_Bd_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D0star_Bd_formula")+name_suffix_vector[i]).c_str(), "@0*(1.-@1)", RooArgList(*f_broad_Bd, *f_D1prime_broad_Bd));  
    f_D0star_Bd_formula_vector.push_back(f_D0star_Bd_formula);

    char hNBd2Dststmunu_D1_name[256];
    strcpy (hNBd2Dststmunu_D1_name,"hNBd2Dststmunu_D1");
    strcat (hNBd2Dststmunu_D1_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bd2Dstst_D1;
    if (cfg->UseHAMMER()) {
      Bd2Dstst_D1 = new RooStats::HistFactory::Sample(hNBd2Dststmunu_D1_name,"HAMMER_inithisto_80_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bd2Dstst_D1_Errors;
      Bd2Dstst_D1_Errors.filename = "hammerinithistos.root";
      Bd2Dstst_D1_Errors.histogramname = "HAMMER_inithisto_80_cat"+isocats[i];

      error_histos[i].emplace(hNBd2Dststmunu_D1_name,Bd2Dstst_D1_Errors); 
    }
    else {
      TH1D* hBd2Dstst_D1 = (TH1D*)fin->Get(hNBd2Dststmunu_D1_name);
      Bd2Dstst_D1 = new RooStats::HistFactory::Sample (hNBd2Dststmunu_D1_name,hNBd2Dststmunu_D1_name,datafile_MC);
      Bd2Dstst_D1->AddNormFactor((std::string("MCBd2Dstst_D1")+name_suffix_vector[i]).c_str(),1.0/hBd2Dstst_D1->Integral(),0.0,1.0);
    
      HistogramErrors Bd2Dstst_D1_Errors;
      Bd2Dstst_D1_Errors.filename = cfg->templates();
      Bd2Dstst_D1_Errors.histogramname = hNBd2Dststmunu_D1_name;

      error_histos[i].emplace(hNBd2Dststmunu_D1_name,Bd2Dstst_D1_Errors);
    }
    Bd2Dstst_D1->SetNormalizeByTheory(kFALSE);
    Bd2Dstst_D1->AddNormFactor(toyparams_vector[i]->NBd2Dstst.get_name(),toyparams_vector[i]->NBd2Dstst.get_value(),toyparams_vector[i]->NBd2Dstst.get_min(),toyparams_vector[i]->NBd2Dstst.get_max());
    Bd2Dstst_D1->AddNormFactor((std::string("f_narrow_D1_Bd")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bd2Dstst_D1->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bd_D1.get_name(), toyparams_vector[i]->TF_Dstst_Bd_D1.get_value(), toyparams_vector[i]->TF_Dstst_Bd_D1.get_min(), toyparams_vector[i]->TF_Dstst_Bd_D1.get_max());
    //if (correct_for_sweights) {Bd2Dstst_D1->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D1_Bd_formula = new RooFormulaVar((std::string("f_D1_Bd_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D1_Bd_formula")+name_suffix_vector[i]).c_str(), "(1.-@0)*@1", RooArgList(*f_broad_Bd, *f_D1_narrow_Bd));
    f_D1_Bd_formula_vector.push_back(f_D1_Bd_formula);

    char hNBd2Dststmunu_D2star_name[256];
    strcpy (hNBd2Dststmunu_D2star_name,"hNBd2Dststmunu_D2star");
    strcat (hNBd2Dststmunu_D2star_name,name_suffix_vector[i].c_str());
    RooStats::HistFactory::Sample* Bd2Dstst_D2star;
    if (cfg->UseHAMMER()) {
      Bd2Dstst_D2star = new RooStats::HistFactory::Sample(hNBd2Dststmunu_D2star_name,"HAMMER_inithisto_70_cat"+isocats[i],datafile_MC_hammer_init);
    
      HistogramErrors Bd2Dstst_D2star_Errors;
      Bd2Dstst_D2star_Errors.filename = "hammerinithistos.root";
      Bd2Dstst_D2star_Errors.histogramname = "HAMMER_inithisto_70_cat"+isocats[i];

      error_histos[i].emplace(hNBd2Dststmunu_D2star_name,Bd2Dstst_D2star_Errors); 
    }
    else {
      TH1D* hBd2Dstst_D2star = (TH1D*)fin->Get(hNBd2Dststmunu_D2star_name);
      Bd2Dstst_D2star = new RooStats::HistFactory::Sample (hNBd2Dststmunu_D2star_name,hNBd2Dststmunu_D2star_name,datafile_MC);
      Bd2Dstst_D2star->AddNormFactor((std::string("MCBd2Dstst_D2star")+name_suffix_vector[i]).c_str(),1.0/hBd2Dstst_D2star->Integral(),0.0,1.0);
    
      HistogramErrors Bd2Dstst_D2star_Errors;
      Bd2Dstst_D2star_Errors.filename = cfg->templates();
      Bd2Dstst_D2star_Errors.histogramname = hNBd2Dststmunu_D2star_name;

      error_histos[i].emplace(hNBd2Dststmunu_D2star_name,Bd2Dstst_D2star_Errors);
    }
    Bd2Dstst_D2star->SetNormalizeByTheory(kFALSE);
    Bd2Dstst_D2star->AddNormFactor(toyparams_vector[i]->NBd2Dstst.get_name(),toyparams_vector[i]->NBd2Dstst.get_value(),toyparams_vector[i]->NBd2Dstst.get_min(),toyparams_vector[i]->NBd2Dstst.get_max());
    Bd2Dstst_D2star->AddNormFactor((std::string("f_narrow_D2star_Bd")+name_suffix_vector[i]).c_str(), 0.25, 0.,1.);
    Bd2Dstst_D2star->AddNormFactor(toyparams_vector[i]->TF_Dstst_Bd_D2.get_name(), toyparams_vector[i]->TF_Dstst_Bd_D2.get_value(), toyparams_vector[i]->TF_Dstst_Bd_D2.get_min(), toyparams_vector[i]->TF_Dstst_Bd_D2.get_max());
    //if (correct_for_sweights) {Bd2Dstst_D2star->AddNormFactor(sweight_normfactor.at(i));}
    RooFormulaVar* f_D2star_Bd_formula = new RooFormulaVar((std::string("f_D2star_Bd_formula")+name_suffix_vector[i]).c_str(), (std::string("f_D2star_Bd_formula")+name_suffix_vector[i]).c_str(), "(1.-@0)*(1.-@1)", RooArgList(*f_broad_Bd, *f_D1_narrow_Bd));
    f_D2star_Bd_formula_vector.push_back(f_D2star_Bd_formula);


    //------------------------B02DststHigher--------------------------
    char hNBd2DststmunuHigher_name[256];
    strcpy (hNBd2DststmunuHigher_name,"hNBd2DststmunuHigher");
    strcat (hNBd2DststmunuHigher_name,name_suffix_vector[i].c_str());
    TH1D* hBd2DststHigher = (TH1D*)fin->Get(hNBd2DststmunuHigher_name);
    RooStats::HistFactory::Sample* Bd2DststHigher = new RooStats::HistFactory::Sample(hNBd2DststmunuHigher_name,hNBd2DststmunuHigher_name,datafile_MC);
    Bd2DststHigher->SetNormalizeByTheory(kFALSE);
    Bd2DststHigher->AddNormFactor((std::string("MCBd2DststHigher")+name_suffix_vector[i]).c_str(),1.0/hBd2DststHigher->Integral(),0.0,1.0);
    Bd2DststHigher->AddNormFactor(toyparams_vector[i]->NBd2DststHigher.get_name(),toyparams_vector[i]->NBd2DststHigher.get_value(),toyparams_vector[i]->NBd2DststHigher.get_min(),toyparams_vector[i]->NBd2DststHigher.get_max());
    //if (correct_for_sweights) {Bd2DststHigher->AddNormFactor(sweight_normfactor.at(i));}
    if(cfg->varyshapes())
    {
      char hNBd2DststmunuHigher_m1sigma_name[256];
      strcpy (hNBd2DststmunuHigher_m1sigma_name,"hNBd2DststmunuHigher_m1sigma");
      strcat (hNBd2DststmunuHigher_m1sigma_name,name_suffix_vector[i].c_str());
      char hNBd2DststmunuHigher_p1sigma_name[256];
      strcpy (hNBd2DststmunuHigher_p1sigma_name,"hNBd2DststmunuHigher_p1sigma");
      strcat (hNBd2DststmunuHigher_p1sigma_name,name_suffix_vector[i].c_str());
      Bd2DststHigher->AddHistoSys(std::string("hNBd2DststmunuHigher_variation"),hNBd2DststmunuHigher_m1sigma_name,cfg->templates(),"",hNBd2DststmunuHigher_p1sigma_name,cfg->templates(),"");
    }

    HistogramErrors Bd2DststHigher_Errors;
    Bd2DststHigher_Errors.filename = cfg->templates();
    Bd2DststHigher_Errors.histogramname = hNBd2DststmunuHigher_name;

    error_histos[i].emplace(hNBd2DststmunuHigher_name, Bd2DststHigher_Errors);
    


    //----------------------------LbDLc-------------------------------
    char hNLbDLc_name[256];
    strcpy (hNLbDLc_name,"hNLbDLc");
    strcat (hNLbDLc_name,name_suffix_vector[i].c_str());
    TH1D* hLbDLc = (TH1D*)fin->Get(hNLbDLc_name);
    RooStats::HistFactory::Sample* LbDLc = new RooStats::HistFactory::Sample(hNLbDLc_name,hNLbDLc_name,datafile_MC);
    LbDLc->SetNormalizeByTheory(kFALSE);
    LbDLc->AddNormFactor((std::string("MCLbDLc")+name_suffix_vector[i]).c_str(),1.0/hLbDLc->Integral(),0.0,1.0);
    LbDLc->AddNormFactor(toyparams_vector[i]->NLbDLc.get_name(),toyparams_vector[i]->NLbDLc.get_value(),toyparams_vector[i]->NLbDLc.get_min(),toyparams_vector[i]->NLbDLc.get_max());
    //if (correct_for_sweights) {LbDLc->AddNormFactor(sweight_normfactor.at(i));}
    
    HistogramErrors LbDLc_Errors;
    LbDLc_Errors.filename = cfg->templates();
    LbDLc_Errors.histogramname = hNLbDLc_name;

    error_histos[i].emplace(hNLbDLc_name, LbDLc_Errors);
    
    //-------------------MisID----------------------------------------
    char hNMisID_name[256];
    strcpy (hNMisID_name,"hNMisID");
    strcat (hNMisID_name,name_suffix_vector[i].c_str());
    TH1D* hMisID = (TH1D*)fin->Get(hNMisID_name);
    RooStats::HistFactory::Sample* MisID = new RooStats::HistFactory::Sample(hNMisID_name,hNMisID_name,datafile_MC);
    MisID->SetNormalizeByTheory(kFALSE);
    MisID->AddNormFactor((std::string("MCMisID")+name_suffix_vector[i]).c_str(),1.0/hMisID->Integral(),0.0,1.0); //Fixed
    MisID->AddNormFactor(toyparams_vector[i]->NMisID.get_name(),toyparams_vector[i]->NMisID.get_value(),toyparams_vector[i]->NMisID.get_min(),toyparams_vector[i]->NMisID.get_max());
    //if (correct_for_sweights) {MisID->AddNormFactor(sweight_normfactor.at(i));}
    if(cfg->ConstrainDataComps()){
      MisID->AddOverallSys(std::string("NormUnc_MisID")+name_suffix_vector[i].c_str(), toyparams_vector[i]->MisID_NormUnc_down.get_value(), toyparams_vector[i]->MisID_NormUnc_up.get_value()); 
      //in this case NMisID is fixed
    }

    HistogramErrors MisID_Errors;
    MisID_Errors.filename = cfg->templates();
    MisID_Errors.histogramname = hNMisID_name;
    
    error_histos[i].emplace(hNMisID_name, MisID_Errors);
  
    //------------------WS--------------------------------------------
    char hNWS_name[256];
    strcpy (hNWS_name,"hNWS");
    strcat (hNWS_name,name_suffix_vector[i].c_str());
    TH1D* hWS = (TH1D*)fin->Get(hNWS_name);
    RooStats::HistFactory::Sample* WS = new RooStats::HistFactory::Sample(hNWS_name,hNWS_name,datafile_MC);
    WS->SetNormalizeByTheory(kFALSE);
    WS->AddNormFactor((std::string("MCWS")+name_suffix_vector[i]).c_str(),1.0/hWS->Integral(),0.0,1.0); //Fixed
    WS->AddNormFactor(toyparams_vector[i]->NWS.get_name(),toyparams_vector[i]->NWS.get_value(),toyparams_vector[i]->NWS.get_min(),toyparams_vector[i]->NWS.get_max());
    //if (correct_for_sweights) {WS->AddNormFactor(sweight_normfactor.at(i));}
    if(cfg->ConstrainDataComps()){
      WS->AddOverallSys(std::string("NormUnc_WS")+name_suffix_vector[i].c_str(), toyparams_vector[i]->WS_NormUnc_down.get_value(), toyparams_vector[i]->WS_NormUnc_up.get_value());
      WS->AddOverallSys(std::string("RS_WS_ratio_WS").c_str(),                     toyparams_vector[i]->RS_WS_ratio_down.get_value(), toyparams_vector[i]->RS_WS_ratio_up.get_value());
      //in this case NWS is fixed
    }

    HistogramErrors WS_Errors;
    WS_Errors.filename = cfg->templates();
    WS_Errors.histogramname = hNWS_name;
    
    error_histos[i].emplace(hNWS_name, WS_Errors);

    if(cfg->MCstat()){
      DplusMuNu->ActivateStatError();
      DstMuNu->ActivateStatError();
      DplusTauNu->ActivateStatError();
      DstTauNu->ActivateStatError();
      Bd2DD_DDs->ActivateStatError();
      Bu2DD_DDs->ActivateStatError();
      BdDDs->ActivateStatError();
      BuDDs->ActivateStatError();
      Bd2DD_MultiBody->ActivateStatError();
      Bu2DD_MultiBody->ActivateStatError();
      Bu2Dstst_D0star->ActivateStatError();
      Bu2Dstst_D2star->ActivateStatError();
      Bu2Dstst_D1->ActivateStatError();
      Bu2Dstst_D1prime->ActivateStatError();
      Bd2Dstst_D0star->ActivateStatError();
      Bd2Dstst_D2star->ActivateStatError();
      Bd2Dstst_D1->ActivateStatError();
      Bd2Dstst_D1prime->ActivateStatError();
      Bd2DststHigher->ActivateStatError();
      LbDLc->ActivateStatError();
      MisID->ActivateStatError();
      WS->ActivateStatError();
    }
    
    std::cout << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    std::cout << "            Adding the samples to channel "                        << isocats[i] << std::endl;
    std::cout << "-------------------------------------------------------------------------------" << std::endl;
    AddSampleToChannel(channels_vector[i], DplusMuNu,       isocats[i], std::string("DplusMuNu"));
    AddSampleToChannel(channels_vector[i], DstMuNu,         isocats[i], std::string("DstMuNu"));
    AddSampleToChannel(channels_vector[i], DplusTauNu,      isocats[i], std::string("DplusTauNu"));
    AddSampleToChannel(channels_vector[i], DstTauNu,        isocats[i], std::string("DstTauNu"));
    AddSampleToChannel(channels_vector[i], Bd2DD_DDs,       isocats[i], std::string("Bd2DD_DDs"));
    AddSampleToChannel(channels_vector[i], Bu2DD_DDs,       isocats[i], std::string("Bu2DD_DDs"));
    AddSampleToChannel(channels_vector[i], BdDDs,           isocats[i], std::string("BdDDs"));
    AddSampleToChannel(channels_vector[i], BuDDs,           isocats[i], std::string("BuDDs"));
    AddSampleToChannel(channels_vector[i], Bd2DD_MultiBody, isocats[i], std::string("Bd2DD_MultiBody"));
    AddSampleToChannel(channels_vector[i], Bu2DD_MultiBody, isocats[i], std::string("Bu2DD_MultiBody"));
    AddSampleToChannel(channels_vector[i], Bu2Dstst_D0star, isocats[i], std::string("Bu2Dstst_D0star"));
    AddSampleToChannel(channels_vector[i], Bu2Dstst_D2star, isocats[i], std::string("Bu2Dstst_D2star"));
    AddSampleToChannel(channels_vector[i], Bu2Dstst_D1,     isocats[i], std::string("Bu2Dstst_D1"));
    AddSampleToChannel(channels_vector[i], Bu2Dstst_D1prime,isocats[i], std::string("Bu2Dstst_D1prime"));
    AddSampleToChannel(channels_vector[i], Bd2Dstst_D0star, isocats[i], std::string("Bd2Dstst_D0star"));
    AddSampleToChannel(channels_vector[i], Bd2Dstst_D2star, isocats[i], std::string("Bd2Dstst_D2star"));
    AddSampleToChannel(channels_vector[i], Bd2Dstst_D1,     isocats[i], std::string("Bd2Dstst_D1"));
    AddSampleToChannel(channels_vector[i], Bd2Dstst_D1prime,isocats[i], std::string("Bd2Dstst_D1prime"));
    AddSampleToChannel(channels_vector[i], Bd2DststHigher,  isocats[i], std::string("Bd2DststHigher"));
    AddSampleToChannel(channels_vector[i], LbDLc,           isocats[i], std::string("LbDLc"));
    AddSampleToChannel(channels_vector[i], MisID,           isocats[i], std::string("MisID"));
    AddSampleToChannel(channels_vector[i], WS,              isocats[i], std::string("WS"));
    std::cout << "-------------------------------------------------------------------------------" << std::endl;

  }
  
  for (unsigned int i=0; i<name_suffix_vector.size(); i++){
    meas.AddChannel(*channels_vector[i]);
  }
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "                            Collecting histograms                              " << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  meas.CollectHistograms();
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "                            Fixing the parameters                              " << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  //Fix the MC normaliztion of the histograms.
  for (unsigned int i=0; i<name_suffix_vector.size(); i++){
    //meas.AddConstantParam(sweight_normfactor.at(i).GetName());
    if (!IsSampleVetoed(isocats[i], std::string("DplusMuNu")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCDplusMu")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCDplusMu")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("DstMuNu")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCDstMu")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCDstMu")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("DplusTauNu")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCDplusTau")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCDplusTau")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("DstTauNu")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCDstTau")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCDstTau")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2DD_DDs"))){
      meas.AddConstantParam((std::string("MCBd2DD_DDs")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2DD_DDs")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2DD_DDs"))){
      meas.AddConstantParam((std::string("MCBu2DD_DDs")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2DD_DDs")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("BdDDs"))){
      meas.AddConstantParam((std::string("MCBdDDs")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBdDDs")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("BuDDs"))){
      meas.AddConstantParam((std::string("MCBuDDs")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBuDDs")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2DD_MultiBody"))){
      meas.AddConstantParam((std::string("MCBd2DD_MultiBody")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2DD_MultiBody")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2DD_MultiBody"))){
      meas.AddConstantParam((std::string("MCBu2DD_MultiBody")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2DD_MultiBody")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D0star")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBu2Dstst_D0star")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2Dstst_D0star")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
        if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D2star")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBu2Dstst_D2star")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2Dstst")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
        if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D1")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBu2Dstst_D1")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2Dstst_D1")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
        if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D1prime")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBu2Dstst_D1prime")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBu2Dstst_D1prime")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D0star")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBd2Dstst_D0star")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2Dstst")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D2star")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBd2Dstst_D2star")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2Dstst")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D1")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBd2Dstst_D1")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2Dstst")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D1prime")) and !cfg->UseHAMMER()){
      meas.AddConstantParam((std::string("MCBd2Dstst_D1prime")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2Dstst")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2DststHigher"))){
      meas.AddConstantParam((std::string("MCBd2DststHigher")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCBd2DststHigher")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("LbDLc"))){
      meas.AddConstantParam((std::string("MCLbDLc")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCLbDLc")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("MisID"))){
      meas.AddConstantParam((std::string("MCMisID")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCMisID")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("WS"))){
      meas.AddConstantParam((std::string("MCWS")+name_suffix_vector[i]).c_str());
      std::cout << "INFO: Parameter " << (std::string("MCWS")+name_suffix_vector[i]).c_str() << " is declared fixed." << std::endl;
    }
  }
  //Fix the raw normalization of the data derived templates if you want it fixed or you want it constrained
  for (unsigned int i=0; i<name_suffix_vector.size(); i++){
    if(!IsSampleVetoed(isocats[i], std::string("MisID"))){
      if(cfg->ConstrainDataComps() || (toyparams_vector[i]->NMisID.get_min() == toyparams_vector[i]->NMisID.get_max())){
        meas.AddConstantParam(toyparams_vector[i]->NMisID.get_name());
        std::cout << "INFO: Parameter "<< std::string(toyparams_vector[i]->NMisID.get_name()) << " is declared fixed." << std::endl;
      }
    }
    if(!IsSampleVetoed(isocats[i], std::string("WS"))){
      if(cfg->ConstrainDataComps() || (toyparams_vector[i]->NWS.get_min() == toyparams_vector[i]->NWS.get_max())){
        meas.AddConstantParam(toyparams_vector[i]->NWS.get_name());
        std::cout << "INFO: Parameter "<< std::string(toyparams_vector[i]->NWS.get_name()) << " is declared fixed." << std::endl;
      }  
    }
  }
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "Parameters fixed" <<std::endl;
  
  //----Change the constraints on some of the systematics
  if(cfg->varyshapes()){
    //TODO: This bit of the code does not work anymore. Understand why
    //meas.AddUniformSyst(std::string("MultiBody_variation"));
    //meas.AddUniformSyst(std::string("MultiBody_quad_variation"));
    //meas.AddNoSyst(std::string("MultiBody_variation")+name_suffix_vector[0]);
    //meas.AddNoSyst(std::string("MultiBody_quad_variation")+name_suffix_vector[0]);
  }
  
  RooWorkspace *w;
  w= RooStats::HistFactory::MakeModelAndMeasurementFast(meas);
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "                            Importing the formulas inside the model            " << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  
  
  std::cout << "Importing "<< RawRDplus_unblind->GetName() << " inside the workspace ";
  std::cout << w->import(*RawRDplus_unblind) << std::endl;
  std::cout << "Importing "<< RawRDst_unblind->GetName() << " inside the workspace "; 
  std::cout << w->import(*RawRDst_unblind) << std::endl;

  for (unsigned int i =0 ; i<name_suffix_vector.size(); i++){
    if (!IsSampleVetoed(isocats[i], std::string("Bd2DD_DDs")) && !IsSampleVetoed(isocats[i], std::string("Bd2DD_MultiBody"))){
      std::cout << "Importing " << fBd2DD_DDs_formula_vector[i]->GetName()<< " inside the workspace "; 
      std::cout << w->import(*fBd2DD_DDs_formula_vector[i]) << std::endl;
      std::cout << "Importing " << fBd2DD_MultiBody_formula_vector[i]->GetName()<< " inside the workspace "; 
      std::cout << w->import(*fBd2DD_MultiBody_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2DD_DDs")) && !IsSampleVetoed(isocats[i], std::string("Bu2DD_MultiBody"))){
      std::cout << "Importing " << fBu2DD_DDs_formula_vector[i]->GetName() << " inside the workspace "; 
      std::cout << w->import(*fBu2DD_DDs_formula_vector[i]) << std::endl;
      std::cout << "Importing " << fBu2DD_MultiBody_formula_vector[i]->GetName() << " inside the workspace "; 
      std::cout << w->import(*fBu2DD_MultiBody_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("BdDDs"))){
      std::cout << "Importing " << fBdDDs_formula_vector[i]->GetName() << " inside the workspace "; 
      std::cout << w->import(*fBdDDs_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("BuDDs"))){
      std::cout << "Importing " << fBuDDs_formula_vector[i]->GetName() << " inside the workspace "; 
      std::cout << w->import(*fBuDDs_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D0star"))){
      std::cout << "Importing " << f_D0star_Bu_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D0star_Bu_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D2star"))){
      std::cout << "Importing " << f_D2star_Bu_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D2star_Bu_formula_vector[i]) << std::endl;
    } 
    if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D1"))){
      std::cout << "Importing " << f_D1_Bu_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D1_Bu_formula_vector[i]) << std::endl;
    } 
    if (!IsSampleVetoed(isocats[i], std::string("Bu2Dstst_D1prime"))){
      std::cout << "Importing " << f_D1prime_Bu_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D1prime_Bu_formula_vector[i]) << std::endl;
    }
        if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D0star"))){
      std::cout << "Importing " << f_D0star_Bd_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D0star_Bd_formula_vector[i]) << std::endl;
    }
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D2star"))){
      std::cout << "Importing " << f_D2star_Bd_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D2star_Bd_formula_vector[i]) << std::endl;
    } 
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D1"))){
      std::cout << "Importing " << f_D1_Bd_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D1_Bd_formula_vector[i]) << std::endl;
    } 
    if (!IsSampleVetoed(isocats[i], std::string("Bd2Dstst_D1prime"))){
      std::cout << "Importing " << f_D1prime_Bd_formula_vector[i]->GetName() << " inside the workspace ";
      std::cout << w->import(*f_D1prime_Bd_formula_vector[i]) << std::endl;
    } 
  }
  
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  RooStats::ModelConfig *mc = (RooStats::ModelConfig*) w->obj("ModelConfig");
  ((RooRealVar*)(mc->GetNuisanceParameters()->find("Lumi")))->setConstant(kTRUE);

  std::vector<RooRealVar*> x_vector;
  std::vector<RooRealVar*> y_vector;
  std::vector<RooRealVar*> z_vector;

  RooArgSet *obs = (RooArgSet*) mc->GetObservables();
  for (unsigned int i=0; i<name_suffix_vector.size(); i++){
    x_vector.push_back((RooRealVar*) obs->find((std::string("obs_x_RDplus_kinematic")+name_suffix_vector[i]).c_str()));
    y_vector.push_back((RooRealVar*) obs->find((std::string("obs_y_RDplus_kinematic")+name_suffix_vector[i]).c_str()));
    z_vector.push_back((RooRealVar*) obs->find((std::string("obs_z_RDplus_kinematic")+name_suffix_vector[i]).c_str()));
    if (cfg->UseHAMMER()) {
      cat_HAMMER_PDFs_00[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_01[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_10[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_11[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_20[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_30[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_40[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_50[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_60[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_70[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_80[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);
      cat_HAMMER_PDFs_90[i].SetObservables(*x_vector[i],*y_vector[i],*z_vector[i]);

    }
  }

  RooCategory *idx = (RooCategory*) obs->find("channelCat");
  RooAbsData *data = (RooAbsData*) w->data("obsData");

  int num_channels = meas.GetChannels().size();
  
  std::cout << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  std::cout << "                            Modifying the model                                " << std::endl;
  std::cout << "-------------------------------------------------------------------------------" << std::endl;

  RooRealVar* RawRDplus_dummy = (RooRealVar*) mc->GetNuisanceParameters()->find("RawRDplus_dummy");
  RooRealVar* RawRDst_dummy   = (RooRealVar*) mc->GetNuisanceParameters()->find("RawRDst_dummy");
  std::vector<RooRealVar*> fBd2DD_DDs_vector;
  std::vector<RooRealVar*> fBu2DD_DDs_vector;
  std::vector<RooRealVar*> fBdDDs_vector;
  std::vector<RooRealVar*> fBuDDs_vector;
  std::vector<RooRealVar*> fBd2DD_MultiBody_vector;
  std::vector<RooRealVar*> fBu2DD_MultiBody_vector;
  std::vector<RooHistFunc*> init_model_Hc0l0_vector;
  std::vector<RooHistFunc*> init_model_Hc0l1_vector;
  std::vector<RooHistFunc*> init_model_Hc1l0_vector;
  std::vector<RooHistFunc*> init_model_Hc1l1_vector;
  std::vector<RooRealVar*> f_broad_D1prime_Bu_vector;
  std::vector<RooRealVar*> f_broad_D0star_Bu_vector;
  std::vector<RooRealVar*> f_narrow_D1_Bu_vector;
  std::vector<RooRealVar*> f_narrow_D2star_Bu_vector;
  std::vector<RooRealVar*> f_broad_D1prime_Bd_vector;
  std::vector<RooRealVar*> f_broad_D0star_Bd_vector;
  std::vector<RooRealVar*> f_narrow_D1_Bd_vector;
  std::vector<RooRealVar*> f_narrow_D2star_Bd_vector;
  
  
  for (unsigned int i=0; i<num_channels; i++){    
    //The parameters are pushed even though they are empty. In those cases they will be Null pointers
    fBd2DD_DDs_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBd2DD_DDs")+name_suffix_vector[i]).c_str()));
    fBu2DD_DDs_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBu2DD_DDs")+name_suffix_vector[i]).c_str()));
    fBdDDs_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBdDDs")+name_suffix_vector[i]).c_str()));
    fBuDDs_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBuDDs")+name_suffix_vector[i]).c_str()));
    fBd2DD_MultiBody_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBd2DD_MultiBody")+name_suffix_vector[i]).c_str()));
    fBu2DD_MultiBody_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("fBu2DD_MultiBody")+name_suffix_vector[i]).c_str()));
    init_model_Hc0l0_vector.push_back((RooHistFunc*) w->obj((std::string("hNDpMuNu")+name_suffix_vector[i]+std::string("_RDplus_kinematic")+name_suffix_vector[i]+std::string("_nominal")).c_str()));
    init_model_Hc0l1_vector.push_back((RooHistFunc*) w->obj((std::string("hNDpTauNu")+name_suffix_vector[i]+std::string("_RDplus_kinematic")+name_suffix_vector[i]+std::string("_nominal")).c_str()));
    init_model_Hc1l0_vector.push_back((RooHistFunc*) w->obj((std::string("hNDstMuNu")+name_suffix_vector[i]+std::string("_RDplus_kinematic")+name_suffix_vector[i]+std::string("_nominal")).c_str()));
    init_model_Hc1l1_vector.push_back((RooHistFunc*) w->obj((std::string("hNDstTauNu")+name_suffix_vector[i]+std::string("_RDplus_kinematic")+name_suffix_vector[i]+std::string("_nominal")).c_str()));
    f_broad_D1prime_Bu_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_broad_D1prime_Bu")+name_suffix_vector[i]).c_str()));
    f_broad_D0star_Bu_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_broad_D0star_Bu")+name_suffix_vector[i]).c_str()));
    f_narrow_D1_Bu_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_narrow_D1_Bu")+name_suffix_vector[i]).c_str()));
    f_narrow_D2star_Bu_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_narrow_D2star_Bu")+name_suffix_vector[i]).c_str()));
    f_broad_D1prime_Bd_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_broad_D1prime_Bd")+name_suffix_vector[i]).c_str()));
    f_broad_D0star_Bd_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_broad_D0star_Bd")+name_suffix_vector[i]).c_str()));
    f_narrow_D1_Bd_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_narrow_D1_Bd")+name_suffix_vector[i]).c_str()));
    f_narrow_D2star_Bd_vector.push_back((RooRealVar*) mc->GetNuisanceParameters()->find((std::string("f_narrow_D2star_Bd")+name_suffix_vector[i]).c_str()));
  }

  RooSimultaneous *model_ = (RooSimultaneous*)mc->GetPdf(); // Original model
  RooSimultaneous *model  = new RooSimultaneous("simPdf_modified","simPdf_modified", *idx); // New model
  //RooSimultaneous *model_sweightcorr = new RooSimultaneous("simPdf_modified_sweightcorr","simPdf_modified_sweightcorr", *idx); // New model + sWeights correction
  RooCustomizer* cust_;
  RooAbsPdf* pdf_;
  for (unsigned int j = 0; j < num_channels; j++) {
    idx->setIndex(j);
    pdf_ = model_->getPdf(idx->getLabel());
    cust_ = new RooCustomizer(*pdf_,"cust_");
    if (RawRDplus_dummy){
      cust_->replaceArg(*RawRDplus_dummy,*RawRDplus_unblind);
      std::cout << "INFO: Parameter " << RawRDplus_dummy->GetName() << " has been substituted with parameter " << RawRDplus_unblind->GetName() << std::endl; 
    }
    if (RawRDst_dummy){
      cust_->replaceArg(*RawRDst_dummy,*RawRDst_unblind);
      std::cout << "INFO: Parameter " << RawRDst_dummy->GetName() << " has been substituted with parameter " << RawRDst_unblind->GetName() << std::endl; 
    }
    if (cfg->UseHAMMER()) {
      if (init_model_Hc0l0_vector[j]) {
        cust_->replaceArg(*init_model_Hc0l0_vector[j], cat_HAMMER_PDFs_00[j]);
        std::cout << "INFO: Function " << init_model_Hc0l0_vector[j]->GetName() << " has been substituted with function " << cat_HAMMER_PDFs_00[j].GetName() << std::endl;
      }   
      if (init_model_Hc0l1_vector[j]) {
        cust_->replaceArg(*init_model_Hc0l1_vector[j], cat_HAMMER_PDFs_01[j]);
        std::cout << "INFO: Function " << init_model_Hc0l1_vector[j]->GetName() << " has been substituted with function " << cat_HAMMER_PDFs_01[j].GetName() << std::endl;
      }   
      if (init_model_Hc1l0_vector[j]) {
        cust_->replaceArg(*init_model_Hc1l0_vector[j], cat_HAMMER_PDFs_10[j]);
        std::cout << "INFO: Function " << init_model_Hc1l0_vector[j]->GetName() << " has been substituted with function " << cat_HAMMER_PDFs_10[j].GetName() << std::endl;
      }   
      if (init_model_Hc1l1_vector[j]) {
        cust_->replaceArg(*init_model_Hc1l1_vector[j], cat_HAMMER_PDFs_11[j]);
        std::cout << "INFO: Function " << init_model_Hc1l1_vector[j]->GetName() << " has been substituted with function " << cat_HAMMER_PDFs_11[j].GetName() << std::endl;
      }   
    }
    if (fBd2DD_DDs_vector[j]){
      cust_->replaceArg(*fBd2DD_DDs_vector[j],*fBd2DD_DDs_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBd2DD_DDs_vector[j]->GetName() << " has been substituted with parameter " << fBd2DD_DDs_formula_vector[j]->GetName() << std::endl; 
    }
    if (fBu2DD_DDs_vector[j]){
      cust_->replaceArg(*fBu2DD_DDs_vector[j],*fBu2DD_DDs_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBu2DD_DDs_vector[j]->GetName() << " has been substituted with parameter " << fBu2DD_DDs_formula_vector[j]->GetName() << std::endl; 
    }
    if (fBdDDs_vector[j]){
      cust_->replaceArg(*fBdDDs_vector[j],*fBdDDs_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBdDDs_vector[j]->GetName() << " has been substituted with parameter " << fBdDDs_formula_vector[j]->GetName() << std::endl; 
    }
    if (fBuDDs_vector[j]){
      cust_->replaceArg(*fBuDDs_vector[j],*fBuDDs_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBuDDs_vector[j]->GetName() << " has been substituted with parameter " << fBuDDs_formula_vector[j]->GetName() << std::endl; 
    }
    if (fBd2DD_MultiBody_vector[j]){
      cust_->replaceArg(*fBd2DD_MultiBody_vector[j],*fBd2DD_MultiBody_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBd2DD_MultiBody_vector[j]->GetName() << " has been substituted with parameter " << fBd2DD_MultiBody_formula_vector[j]->GetName() << std::endl; 
    }
    if (fBu2DD_MultiBody_vector[j]){
      cust_->replaceArg(*fBu2DD_MultiBody_vector[j],*fBu2DD_MultiBody_formula_vector[j]);
      std::cout << "INFO: Parameter " << fBu2DD_MultiBody_vector[j]->GetName() << " has been substituted with parameter " << fBu2DD_MultiBody_formula_vector[j]->GetName() << std::endl; 
    }
    if (f_broad_D1prime_Bu_vector[j]){
      cust_->replaceArg(*f_broad_D1prime_Bu_vector[j],*f_D1prime_Bu_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_broad_D1prime_Bu_vector[j]->GetName() << " has been substituted with parameter " << f_D1prime_Bu_formula_vector[j]->GetName() << std::endl;
    }
    if (f_broad_D0star_Bu_vector[j]){
      cust_->replaceArg(*f_broad_D0star_Bu_vector[j],*f_D0star_Bu_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_broad_D0star_Bu_vector[j]->GetName() << " has been substituted with parameter " << f_D0star_Bu_formula_vector[j]->GetName() << std::endl;
    }
    if (f_narrow_D1_Bu_vector[j]){
      cust_->replaceArg(*f_narrow_D1_Bu_vector[j],*f_D1_Bu_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_narrow_D1_Bu_vector[j]->GetName() << " has been substituted with parameter " << f_D1_Bu_formula_vector[j]->GetName() << std::endl;
    }
    if (f_narrow_D2star_Bu_vector[j]){
      cust_->replaceArg(*f_narrow_D2star_Bu_vector[j],*f_D2star_Bu_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_narrow_D2star_Bu_vector[j]->GetName() << " has been substituted with parameter " << f_D2star_Bu_formula_vector[j]->GetName() << std::endl;
    }
    if (f_broad_D1prime_Bd_vector[j]){
      cust_->replaceArg(*f_broad_D1prime_Bd_vector[j],*f_D1prime_Bd_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_broad_D1prime_Bd_vector[j]->GetName() << " has been substituted with parameter " << f_D1prime_Bd_formula_vector[j]->GetName() << std::endl;
    }
    if (f_broad_D0star_Bd_vector[j]){
      cust_->replaceArg(*f_broad_D0star_Bd_vector[j],*f_D0star_Bd_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_broad_D0star_Bd_vector[j]->GetName() << " has been substituted with parameter " << f_D0star_Bd_formula_vector[j]->GetName() << std::endl;
    }
    if (f_narrow_D1_Bd_vector[j]){
      cust_->replaceArg(*f_narrow_D1_Bd_vector[j],*f_D1_Bd_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_narrow_D1_Bd_vector[j]->GetName() << " has been substituted with parameter " << f_D1_Bd_formula_vector[j]->GetName() << std::endl;
    }
    if (f_narrow_D2star_Bd_vector[j]){
      cust_->replaceArg(*f_narrow_D2star_Bd_vector[j],*f_D2star_Bd_formula_vector[j]);
      std::cout << "INFO: Parameter " << f_narrow_D2star_Bd_vector[j]->GetName() << " has been substituted with parameter " << f_D2star_Bd_formula_vector[j]->GetName() << std::endl;
    }
    
    model->addPdf((RooAbsPdf&)*cust_->build(),idx->getLabel());

    //if (correct_for_sweights) {
    //  sweightcorr_datahist.push_back(RooDataHist(("sweightcorr_datahist_channel_"+name_suffix_vector[j]).c_str(),("sweightcorr_datahist_channel_"+name_suffix_vector[j]).c_str(),RooArgList(*x_vector[j],*y_vector[j],*z_vector[j]),sweightcorrhisto->at(j)));
    //  sweightcorr_histpdf.push_back(RooHistPdf(("sweightcorr_histpdf_channel_"+name_suffix_vector[j]).c_str(),("sweightcorr_histpdf_channel_"+name_suffix_vector[j]).c_str(),RooArgList(*x_vector[j],*y_vector[j],*z_vector[j]),sweightcorr_datahist.at(j)));
    //  sweightcorr_model_chan.push_back(RooProdPdf(("sweightcorrected_model_chan"+name_suffix_vector[j]).c_str(),("sweightcorrected_model_chan"+name_suffix_vector[j]).c_str(),*model->getPdf(idx->getLabel()),sweightcorr_histpdf.at(j)));
    //  model_sweightcorr->addPdf(sweightcorr_model_chan.at(j),idx->getLabel());
    //}

  }
  
  const RooArgSet* POI_argset_ = mc->GetParametersOfInterest();
  RooArgSet* POI_argset = new RooArgSet(*POI_argset_);
  
  POI_argset->add(*RawRDplus);
  std::cout << "INFO: Parameter " << RawRDplus->GetName() << " has been made POI" << std::endl;
  POI_argset->add(*RawRDst);
  std::cout << "INFO: Parameter " << RawRDst->GetName() << " has been made POI" << std::endl;

  //Define a different starting point for the most correlated parameters
  if (cfg->varyshapes())
  {
    RooRealVar *tempVar; //dummy pointer
    TIterator *iterpars = mc->GetNuisanceParameters()->createIterator();
    while((tempVar = (RooRealVar*) iterpars->Next())) {
      TString theName = tempVar->GetName();
      if(theName.Contains("alpha_MultiBody_variation")){
        tempVar->setVal(0.5);
      }
      if(theName.Contains("alpha_MultiBody_quad_variation")){
        tempVar->setVal(-0.3);
      }
      if(theName.Contains("hNBd2DststmunuHigher_variation")){
        tempVar->setVal(0.1);
      }
    } 
    RooRealVar *tempVar_poi; //dummy pointer
    TIterator *iterpars_pois = mc->GetParametersOfInterest()->createIterator();
    while((tempVar_poi = (RooRealVar*) iterpars_pois->Next())) {
      TString theName = tempVar_poi->GetName();
      if(theName.Contains("alpha_MultiBody_variation")){
        tempVar_poi->setVal(0.5);
      }
      if(theName.Contains("alpha_MultiBody_quad_variation")){
        tempVar_poi->setVal(-0.3);
      }
    } 
  }
  /*
  if (cfg->ConstrainDataComps())
  {
    RooRealVar *tempVar; //dummy pointer
    TIterator *iterpars = mc->GetNuisanceParameters()->createIterator();
    while((tempVar = (RooRealVar*) iterpars->Next())) {
      TString theName = tempVar->GetName();
      if(theName.Contains("alpha_NormUnc_MisID")){
        tempVar->setVal(0.2);
      }
      if(theName.Contains("alpha_NormUnc_WS")){
        tempVar->setVal(0.5);
      }
      if(theName.Contains("alpha_RS_WS_ratio_WS")){
        tempVar->setVal(0.3);
      }
    }
  }
  */



  for (unsigned int i=0; i<num_channels; i++){

          
    if (f_Bu_vector[i]){
      if (toyparams_vector[i]->f_Bu.get_min()==toyparams_vector[i]->f_Bu.get_max()) {f_Bu_vector[i]->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_Bu_vector[i]->GetName() << " is declared constant" << std::endl;
    } 
    if (f_DD_Bu_vector[i]){
      if (toyparams_vector[i]->f_DD_Bu.get_min()==toyparams_vector[i]->f_DD_Bu.get_max()) {f_DD_Bu_vector[i]->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_DD_Bu_vector[i]->GetName() << " is declared constant" << std::endl;
    } 
    if (f_DD_Bd_vector[i]){
      if (toyparams_vector[i]->f_DD_Bd.get_min()==toyparams_vector[i]->f_DD_Bd.get_max()) {f_DD_Bd_vector[i]->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_DD_Bd_vector[i]->GetName() << " is declared constant" << std::endl;
    } 
    if (f_DD_tauonic_Bu_vector[i]){
      if (toyparams_vector[i]->f_DD_tauonic_Bu.get_min()==toyparams_vector[i]->f_DD_tauonic_Bu.get_max()) {f_DD_tauonic_Bu_vector[i]->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_DD_tauonic_Bu_vector[i]->GetName() << " is declared constant" << std::endl;
    } 
    if (f_DD_tauonic_Bd_vector[i]){
      if (toyparams_vector[i]->f_DD_tauonic_Bd.get_min()==toyparams_vector[i]->f_DD_tauonic_Bd.get_max()) {f_DD_tauonic_Bd_vector[i]->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_DD_tauonic_Bd_vector[i]->GetName() << " is declared constant" << std::endl;
    } 
    if (f_broad_Bu){
      if (common_parameters->f_broad_Bu.get_min()==common_parameters->f_broad_Bu.get_max()) {f_broad_Bu->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_broad_Bu->GetName() << " is declared constant" << std::endl;
    }
    if (f_broad_Bd){
      if (common_parameters->f_broad_Bd.get_min()==common_parameters->f_broad_Bd.get_max()) {f_broad_Bd->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_broad_Bd->GetName() << " is declared constant" << std::endl;
    }
    if (f_D1_narrow_Bu){
      if (common_parameters->f_D1_narrow_Bu.get_min()==common_parameters->f_D1_narrow_Bu.get_max()) {f_D1_narrow_Bu->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_D1_narrow_Bu->GetName() << " is declared constant" << std::endl;
    }
    if (f_D1_narrow_Bd){
      if (common_parameters->f_D1_narrow_Bd.get_min()==common_parameters->f_D1_narrow_Bd.get_max()) {f_D1_narrow_Bd->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_D1_narrow_Bd->GetName() << " is declared constant" << std::endl;
    }
    if (f_D1prime_broad_Bu){
      if (common_parameters->f_D1prime_broad_Bu.get_min()==common_parameters->f_D1prime_broad_Bu.get_max()) {f_D1prime_broad_Bu->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_D1prime_broad_Bu->GetName() << " is declared constant" << std::endl;
    }
    if (f_D1prime_broad_Bd){
      if (common_parameters->f_D1prime_broad_Bd.get_min()==common_parameters->f_D1prime_broad_Bd.get_max()) {f_D1prime_broad_Bd->setConstant(kTRUE);}
      std::cout << "INFO: Parameter "<< f_D1prime_broad_Bd->GetName() << " is declared constant" << std::endl;
    }
    if (f_Bu_vector[i]){
      POI_argset->add(*f_Bu_vector[i]);
      std::cout << "INFO: Parameter " << f_Bu_vector[i]->GetName() << " has been made POI" << std::endl; 
    }    
    if (f_DD_Bd_vector[i]){
      POI_argset->add(*f_DD_Bd_vector[i]);
      std::cout << "INFO: Parameter " << f_DD_Bd_vector[i]->GetName() << " has been made POI" << std::endl; 
    }    
    if (f_DD_Bu_vector[i]){
      POI_argset->add(*f_DD_Bu_vector[i]);
      std::cout << "INFO: Parameter " << f_DD_Bu_vector[i]->GetName() << " has been made POI" << std::endl; 
    }    
    if (f_DD_tauonic_Bd_vector[i]){
      POI_argset->add(*f_DD_tauonic_Bd_vector[i]);
      std::cout << "INFO: Parameter " << f_DD_tauonic_Bd_vector[i]->GetName() << " has been made POI" << std::endl; 
    }     
    if (f_DD_tauonic_Bu_vector[i]){
      POI_argset->add(*f_DD_tauonic_Bu_vector[i]);
      std::cout << "INFO: Parameter " << f_DD_tauonic_Bu_vector[i]->GetName() << " has been made POI" << std::endl; 
    }
    if (f_broad_Bu){
      POI_argset->add(*f_broad_Bu);
      std::cout << "INFO: Parameter "<< f_broad_Bu->GetName() << " has been made POI" << std::endl;
    }
    if (f_broad_Bd){
      POI_argset->add(*f_broad_Bd);
      std::cout << "INFO: Parameter "<< f_broad_Bd->GetName() << " has been made POI" << std::endl;
    }
    if (f_D1_narrow_Bu){
      POI_argset->add(*f_D1_narrow_Bu);
      std::cout << "INFO: Parameter "<< f_D1_narrow_Bu->GetName() << " has been made POI" << std::endl;
    }
    if (f_D1_narrow_Bd){
      POI_argset->add(*f_D1_narrow_Bd);
      std::cout << "INFO: Parameter "<< f_D1_narrow_Bd->GetName() << " has been made POI" << std::endl;
    }
    if (f_D1prime_broad_Bu){
      POI_argset->add(*f_D1prime_broad_Bu);
      std::cout << "INFO: Parameter "<< f_D1prime_broad_Bu->GetName() << " has been made POI" << std::endl;
    }
    if (f_D1prime_broad_Bd){
      POI_argset->add(*f_D1prime_broad_Bd);
      std::cout << "INFO: Parameter "<< f_D1prime_broad_Bd->GetName() << " has been made POI" << std::endl;
    }
    if (cfg->UseHAMMER()) {
      if (cfg->ConstrainBtoDBGLParams()) {
        POI_argset->add(*delta_ap0);
        std::cout << "INFO: Parameter "<< delta_ap0->GetName() << " has been made POI" << std::endl; 
      }
      POI_argset->add(*delta_ap1);
      std::cout << "INFO: Parameter "<< delta_ap1->GetName() << " has been made POI" << std::endl; 
      if (cfg->ConstrainBtoDBGLParams()) {
        POI_argset->add(*delta_ap2);
        std::cout << "INFO: Parameter "<< delta_ap2->GetName() << " has been made POI" << std::endl;
      }
      //POI_argset->add(*delta_ap3);
      //std::cout << "INFO: Parameter "<< delta_ap3->GetName() << " has been made POI" << std::endl; 
      if (!cfg->ConstrainBtoDBGLParams()) {
        POI_argset->add(*delta_a00);
        std::cout << "INFO: Parameter "<< delta_a00->GetName() << " has been made POI" << std::endl;
      }
      POI_argset->add(*delta_a01);
      std::cout << "INFO: Parameter "<< delta_a01->GetName() << " has been made POI" << std::endl; 
      if (cfg->ConstrainBtoDBGLParams()) {
        POI_argset->add(*delta_a02);
        std::cout << "INFO: Parameter "<< delta_a02->GetName() << " has been made POI" << std::endl; 
      }
      //POI_argset->add(*delta_a03);
      //std::cout << "INFO: Parameter "<< delta_a03->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_a0);
      std::cout << "INFO: Parameter "<< delta_a0->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_a1);
      std::cout << "INFO: Parameter "<< delta_a1->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_a2);
      std::cout << "INFO: Parameter "<< delta_a2->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_b0);
      std::cout << "INFO: Parameter "<< delta_b0->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_b1);
      std::cout << "INFO: Parameter "<< delta_b1->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_b2);
      std::cout << "INFO: Parameter "<< delta_b2->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_c1);
      std::cout << "INFO: Parameter "<< delta_c1->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_c2);
      std::cout << "INFO: Parameter "<< delta_c2->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_d0);
      std::cout << "INFO: Parameter "<< delta_d0->GetName() << " has been made POI" << std::endl; 
      POI_argset->add(*delta_d1);
      std::cout << "INFO: Parameter "<< delta_d1->GetName() << " has been made POI" << std::endl; 
    }
    
  }
  std::cout << "-------------------------------------------------------------------------------" << std::endl;
  mc->SetParametersOfInterest(*POI_argset);
  
  std::cout << "Model modified" << std::endl;
  
  RooStats::HistFactory::HistFactorySimultaneous* model_hf;
  //RooStats::HistFactory::HistFactorySimultaneous* model_hf_sweightcorr;
  
  //Make sure you are setting the initial value of the model from what you are given in the 
  LoadStartValues(mc, common_parameters, toyparams_vector);

  /*// Creation of the datasets for plotting (with the correct bin errors)
  RooRealVar* bin_weight_plotting = new RooRealVar("bin_weight_plotting","bin_weight_plotting",1.);
  double xbin_center_plotting, ybin_center_plotting, zbin_center_plotting, plotting_bin_content, plotting_bin_error;
  std::vector<RooDataSet> plotting_dataset_vector;
  RooArgSet obs_and_weight_plotting = RooArgSet(*bin_weight_plotting);
  for (unsigned int i=0; i<num_channels; i++){
    obs_and_weight_plotting.add(*x_vector[i]);
    obs_and_weight_plotting.add(*y_vector[i]);
    obs_and_weight_plotting.add(*z_vector[i]);
    plotting_dataset_vector.push_back(RooDataSet(("plotting_dataset_channel_"+std::to_string(i)).c_str(),("plotting_dataset_channel_"+std::to_string(i)).c_str(),RooArgSet(*x_vector[i],*y_vector[i],*z_vector[i],*bin_weight_plotting),RooFit::WeightVar(*bin_weight_plotting),RooFit::StoreError(RooArgSet(*bin_weight_plotting))));
    TH3D* datahisto_;
    TFile* datatemplatefile = new TFile(cfg->templates(),"read");
    datahisto_ = (TH3D*) datatemplatefile->Get(("hNData"+name_suffix_vector[i]).c_str());
    datahisto_->SetDirectory(0);
    datatemplatefile->Close();
    for (unsigned int xbin_idx=0; xbin_idx<x_vector[0]->getBinning().numBins(); xbin_idx++){
      xbin_center_plotting = x_vector[0]->getBinning().binCenter(xbin_idx);
      for (unsigned int ybin_idx=0; ybin_idx<y_vector[0]->getBinning().numBins(); ybin_idx++){
        ybin_center_plotting = y_vector[0]->getBinning().binCenter(ybin_idx);
        for (unsigned int zbin_idx=0; zbin_idx<z_vector[0]->getBinning().numBins(); zbin_idx++){
          zbin_center_plotting = z_vector[0]->getBinning().binCenter(zbin_idx);
          x_vector[i]->setVal(xbin_center_plotting);
          y_vector[i]->setVal(ybin_center_plotting);
          z_vector[i]->setVal(zbin_center_plotting);
          plotting_bin_content = datahisto_->GetBinContent(datahisto_->FindBin(xbin_center_plotting,ybin_center_plotting,zbin_center_plotting));
          plotting_bin_error = datahisto_->GetBinError(datahisto_->FindBin(xbin_center_plotting,ybin_center_plotting,zbin_center_plotting));
          bin_weight_plotting->setVal(plotting_bin_content);
          plotting_dataset_vector[i].add(RooArgSet(*x_vector[i],*y_vector[i],*z_vector[i],*bin_weight_plotting),plotting_bin_content,plotting_bin_error);
        }
      }
    }
  }

  std::map<std::string,RooDataSet*> dataset_map_plotting;
  for (unsigned int i=0; i<num_channels; i++){
    idx->setIndex(i);
    dataset_map_plotting.insert(std::pair<std::string,RooDataSet*>(std::string(idx->getLabel()),&plotting_dataset_vector[i]));
  }

  RooDataSet* plotting_dataset_combined = new RooDataSet("plotting_dataset_combined","plotting_dataset_combined",obs_and_weight_plotting,RooFit::Index(*idx),RooFit::Import(dataset_map_plotting),RooFit::WeightVar(*bin_weight_plotting),RooFit::StoreError(RooArgSet(*bin_weight_plotting)));*/

  // Creation of the data datahists for plotting, and the inverse of the sWeight correction for the model.
  std::vector<TH3D> plotting_TH3D_vector;
  std::vector<RooDataHist> plotting_datahist_vector;
  //std::vector<RooDataHist> sweightcorr_inv_datahist;
  //std::vector<RooHistPdf> sweightcorr_inv_histpdf;
  for (unsigned int i=0; i<num_channels; i++){
    TFile* datatemplatefile = new TFile(cfg->templates(),"read");
    plotting_TH3D_vector.push_back(*((TH3D*) (datatemplatefile->Get(("hNData"+name_suffix_vector[i]).c_str()))->Clone(("TH1"+name_suffix_vector[i]).c_str())));
    datatemplatefile->Close();
    plotting_datahist_vector.push_back(RooDataHist(("plotting_datahist_channel"+name_suffix_vector[i]).c_str(),("plotting_datahist_channel"+name_suffix_vector[i]).c_str(),RooArgList(*x_vector[i],*y_vector[i],*z_vector[i]),&plotting_TH3D_vector[i]));
    //if (correct_for_sweights) {
    //  sweightcorr_inv_datahist.push_back(RooDataHist(("sweightcorr_inv_datahist_channel"+name_suffix_vector[i]).c_str(),("sweightcorr_inv_datahist_channel"+name_suffix_vector[i]).c_str(),RooArgList(*x_vector[i],*y_vector[i],*z_vector[i]),sweightcorrhisto_inv->at(i)));
    //  sweightcorr_inv_histpdf.push_back(RooHistPdf(("sweightcorr_inv_histpdf_channel"+name_suffix_vector[i]).c_str(),("sweightcorr_inv_histpdf_channel"+name_suffix_vector[i]).c_str(),RooArgList(*x_vector[i],*y_vector[i],*z_vector[i]),sweightcorr_inv_datahist.at(i)));
    //}
  }
  std::map<std::string,RooDataHist*> datahist_map_plotting;
  for (unsigned int i=0; i<num_channels; i++){
    idx->setIndex(i);
    datahist_map_plotting.insert(std::pair<std::string,RooDataHist*>(std::string(idx->getLabel()),&plotting_datahist_vector[i]));
  }

  // Preparation of a toy dataset, part 1
  RooRealVar* bin_weight = new RooRealVar("bin_weight","bin_weight",1.);
  std::vector<RooDataSet> toy_dataset_vector;
  for (unsigned int i=0; i<num_channels; i++){
    toy_dataset_vector.push_back(RooDataSet(("toy_dataset_channel"+name_suffix_vector[i]).c_str(),("toy_dataset_channel"+name_suffix_vector[i]).c_str(),RooArgSet(*x_vector[i],*y_vector[i],*z_vector[i],*bin_weight),RooFit::WeightVar(*bin_weight),RooFit::StoreError(RooArgSet(*bin_weight))));
  }
  RooArgSet obs_and_weight = RooArgSet(*bin_weight);
  RooDataSet* toy_dataset_combined;
  std::map<std::string,RooDataSet*> dataset_map;

  if(cfg->justplot()==0) {
    
    model_hf = new RooStats::HistFactory::HistFactorySimultaneous(*model);
    //if (correct_for_sweights) {model_hf_sweightcorr = new RooStats::HistFactory::HistFactorySimultaneous(*model_sweightcorr);}

    std::cout << "Model successfully compiled" << std::endl;
    
  // Preparation of a toy dataset, part 2
    if (cfg->ToyFit()) {
      TRandom randomizer(0);
      double xbin_center, ybin_center, zbin_center, bin_yield, random_bin_yield, bin_yield_error, channel_yield, pdf_integral;
      for (unsigned int i=0; i<num_channels; i++){
        idx->setIndex(i);
        obs_and_weight.add(*x_vector[i]);
        obs_and_weight.add(*y_vector[i]);
        obs_and_weight.add(*z_vector[i]);
        channel_yield = (((RooDataSet*) (w->obj("asimovData")))->reduce(TString("channelCat==channelCat::")+TString(idx->getLabel())))->sumEntries();
        pdf_integral = 0.;
        // Computing the integral of the PDF
        for (unsigned int xbin_idx=0; xbin_idx<x_vector[0]->getBinning().numBins(); xbin_idx++){
          xbin_center = x_vector[0]->getBinning().binCenter(xbin_idx);
          for (unsigned int ybin_idx=0; ybin_idx<y_vector[0]->getBinning().numBins(); ybin_idx++){
            ybin_center = y_vector[0]->getBinning().binCenter(ybin_idx);
            for (unsigned int zbin_idx=0; zbin_idx<z_vector[0]->getBinning().numBins(); zbin_idx++){
              zbin_center = z_vector[0]->getBinning().binCenter(zbin_idx);
              x_vector[i]->setVal(xbin_center);
              y_vector[i]->setVal(ybin_center);
              z_vector[i]->setVal(zbin_center);
              pdf_integral += ((RooAbsPdf*) model_->getPdf(idx->getLabel()))->getVal();
            }
          }
        }
        // Computing the yields per bin
        for (unsigned int xbin_idx=0; xbin_idx<x_vector[0]->getBinning().numBins(); xbin_idx++){
          xbin_center = x_vector[0]->getBinning().binCenter(xbin_idx);
          for (unsigned int ybin_idx=0; ybin_idx<y_vector[0]->getBinning().numBins(); ybin_idx++){
            ybin_center = y_vector[0]->getBinning().binCenter(ybin_idx);
            for (unsigned int zbin_idx=0; zbin_idx<z_vector[0]->getBinning().numBins(); zbin_idx++){
              zbin_center = z_vector[0]->getBinning().binCenter(zbin_idx);
              x_vector[i]->setVal(xbin_center);
              y_vector[i]->setVal(ybin_center);
              z_vector[i]->setVal(zbin_center);
              bin_yield = channel_yield*((RooAbsPdf*) model_hf->getPdf(idx->getLabel()))->getVal()/pdf_integral;
              if (bin_yield>=1) {
                  random_bin_yield = randomizer.Poisson(bin_yield);
                  bin_yield_error = sqrt(bin_yield);
              }
              else {
                  random_bin_yield = 0.;
                  bin_yield_error = 0.;
              }
              bin_weight->setVal(random_bin_yield);
              toy_dataset_vector[i].add(RooArgSet(*x_vector[i],*y_vector[i],*z_vector[i],*bin_weight),random_bin_yield,bin_yield_error);
            }
          }
        }
        dataset_map.insert(std::pair<std::string,RooDataSet*>(std::string(idx->getLabel()),&toy_dataset_vector[i]) );
      }
      toy_dataset_combined = new RooDataSet("toy_dataset_combined","toy_dataset_combined",obs_and_weight,RooFit::Index(*idx),RooFit::Import(dataset_map),RooFit::WeightVar(*bin_weight),RooFit::StoreError(RooArgSet(*bin_weight)));
    }

    RooArgSet set_of_constraints = RooArgSet();
    if (cfg->UseHAMMER()) {
        double nsigma_coverage = 5.; // To be tuned if needed
        if (cfg->FloatFFcoefs_BtoD() and cfg->ConstrainBtoDBGLParams()) {
            std::cout << "INFO: gaussian constraint applied on the BtoD BGL coefficients." << std::endl;
            delta_ap0->setMin(-nsigma_coverage*BtoDBGL_FFparams_sigma_array[0]);
            delta_ap1->setMin(-nsigma_coverage*BtoDBGL_FFparams_sigma_array[1]);
            delta_ap2->setMin(-nsigma_coverage*BtoDBGL_FFparams_sigma_array[2]);
            delta_a01->setMin(-nsigma_coverage*BtoDBGL_FFparams_sigma_array[3]);
            delta_a02->setMin(-nsigma_coverage*BtoDBGL_FFparams_sigma_array[4]);
            delta_ap0->setMax(nsigma_coverage*BtoDBGL_FFparams_sigma_array[0]);
            delta_ap1->setMax(nsigma_coverage*BtoDBGL_FFparams_sigma_array[1]);
            delta_ap2->setMax(nsigma_coverage*BtoDBGL_FFparams_sigma_array[2]);
            delta_a01->setMax(nsigma_coverage*BtoDBGL_FFparams_sigma_array[3]);
            delta_a02->setMax(nsigma_coverage*BtoDBGL_FFparams_sigma_array[4]);
            set_of_constraints.add(*delta_ap0_constraint);
            set_of_constraints.add(*delta_ap1_constraint);
            set_of_constraints.add(*delta_ap2_constraint);
            set_of_constraints.add(*delta_a01_constraint);
            set_of_constraints.add(*delta_a02_constraint);
        }
        if (cfg->FloatFFcoefs_BtoDst() and cfg->ConstrainBtoDstBGLParams()) {
            std::cout << "INFO: gaussian constraint applied on the BtoDst BGL coefficients." << std::endl;
            delta_a0->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[0]);
            delta_a1->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[1]);
            delta_b0->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[2]);
            delta_b1->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[3]);
            delta_c1->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[4]);
            delta_c2->setMin(-nsigma_coverage*BtoDstBGL_FFparams_sigma_array[5]);
            delta_a0->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[0]);
            delta_a1->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[1]);
            delta_b0->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[2]);
            delta_b1->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[3]);
            delta_c1->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[4]);
            delta_c2->setMax(nsigma_coverage*BtoDstBGL_FFparams_sigma_array[5]);
            set_of_constraints.add(*delta_a0_constraint);
            set_of_constraints.add(*delta_a1_constraint);
            set_of_constraints.add(*delta_b0_constraint);
            set_of_constraints.add(*delta_b1_constraint);
            set_of_constraints.add(*delta_c1_constraint);
            set_of_constraints.add(*delta_c2_constraint);
        }
    }

    RooAbsReal* nll_hf;
    if (cfg->ToyFit()==0){
        //if (correct_for_sweights) {nll_hf = model_hf_sweightcorr->createNLL(*data,RooFit::ExternalConstraints(set_of_constraints),RooFit::Offset(kTRUE),RooFit::NumCPU(cfg->FitNumCPU()));}
        //else {nll_hf = model_hf->createNLL(*data,RooFit::ExternalConstraints(set_of_constraints),RooFit::Offset(kTRUE),RooFit::NumCPU(cfg->FitNumCPU()));}
        nll_hf = model_hf->createNLL(*data,RooFit::ExternalConstraints(set_of_constraints),RooFit::Offset(kTRUE),RooFit::NumCPU(cfg->FitNumCPU()));
    }
    else {nll_hf = model_hf->createNLL(*toy_dataset_combined,RooFit::ExternalConstraints(set_of_constraints),RooFit::Offset(kTRUE),RooFit::NumCPU(cfg->FitNumCPU()));}

    RooMinuit* minuit_hf = new RooMinuit(*nll_hf);

    minuit_hf->setErrorLevel(0.5);
    minuit_hf->setStrategy(2);
    minuit_hf->fit("smh"); // First fit, to set up some fine-tuned starting values for the parameters

    std::clock_t start;
    double duration;
    start = std::clock();

    minuit_hf->fit("smh"); // Nominal fit.

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout << "Fitting time: "<< duration << " s." << std::endl;

    /*
    RooFitResult* tempResult=minuit_hf->save("TempResult", "TempResult");

    TFile* output_simpleFit = new TFile(TString(cfg->outputFolder())+TString("/")+TString(cfg->output()), "RECREATE");
    output_simpleFit->cd();

    //--Save the variables and some related quantities
    for (unsigned int j=0; j<toyparams_vector.size();j++){
      for (unsigned int k=0; k<toyparams_vector.at(j)->nparameters(); k++){
        parameter* param = toyparams_vector.at(j)->get_parameter(k);
        std::string parname(param->get_name());
        std::string pardesc(param->get_description());
        TTree* t = new TTree(parname.c_str(), pardesc.c_str());
        double value, error, error_up, error_down, start_value;
        int migrad, status_cov, toynumber;
        t->Branch("value",&value,"value/D");
        t->Branch("error",&error,"error/D");
        t->Branch("start_value",&start_value,"start_value/D");
        t->Branch("migrad",&migrad,"migrad/I");
        t->Branch("status_cov",&status_cov,"status_cov/I"); 

        RooRealVar* FinalParameter = (RooRealVar*)mc->GetParametersOfInterest()->find(param->get_name().c_str());
        if (param->is_blind()){
          value = FinalParameter->getVal() + param->get_blinding_delta();
        }
        else{
          value = FinalParameter->getVal();
        }
        error = FinalParameter->getError();
        if (param->is_blind()){
          start_value = param->get_start_value() + param->get_blinding_delta();
        }
        else{
          start_value = param->get_start_value();
        }

        migrad     = tempResult->status();  //0 means last step (hess in our case) has converged
        status_cov = tempResult->covQual(); //-1 = not available (inversion failes or Hesse failed)
                                            // 0 = available but not positive defined
                                            // 1 = covariance matrix only approximate
                                            // 2 = full matrix but forced pos def
                                            // 3 = full accurate  matrix

        t->Branch("toynumber", &toynumber, "toynumber/I");
        toynumber=cfg->toyn();
        t->Fill();
        t->Write();
        delete t;
      }
    }

    TTree* t = new TTree("migrad","migrad");
    int migrad, status_cov;
    t->Branch("migrad",&migrad,"migrad/I");
    t->Branch("status_cov",&status_cov,"status_cov/I");
    migrad     = tempResult->status(); //0 means last step (in this case hesse) has converged
    
    status_cov = tempResult->covQual(); //-1 = not available (inversion failed or Hesse failed)
                                        // 0 = available but not positive defined
                                        // 1 = covariance only approximate
                                        // 2 = full matrix but forced pos def
                                        // 3 = full accurate matrix

    t->Fill();
    t->Write();

    //Save the correlation between RD and RDst
    for (unsigned int i=0; i<name_suffix_vector.size(); i++){
      std::string correlation_name = std::string("correlation")+name_suffix_vector[i];
      TTree* t_corr = new TTree(correlation_name.c_str(), correlation_name.c_str());
      double RD_RDst_correlation;
      t_corr->Branch("migrad",              &migrad,      "migrad/I"); 
      t_corr->Branch("status_cov",          &status_cov,  "status_cov/I");
      t_corr->Branch("RD_RDst_correlation", &RD_RDst_correlation,"RD_RDst_correlation/D");

      RD_RDst_correlation = tempResult->correlation((std::string("RawRDplus")+name_suffix_vector[i]).c_str(), (std::string("RawRDst")+name_suffix_vector[i]).c_str());
      t_corr->Fill();
      t_corr->Write(); 
    }

    TTree* t_toy = new TTree("toy_number", "toy_number");
    int toynumber;
    t_toy->Branch("toynumber", &toynumber, "toynumber/I");
    toynumber=cfg->toyn();
    t_toy->Fill();
    t_toy->Write();

    // This part is very expensive! Keeping this commented out for the moment
    //---For the moment I am leaving it commented
    
    /*
    TFile* fnll = new TFile("nllplots.root", "RECREATE");
    for (unsigned int j = 0; j <toyparams_vector.size(); j++) {
        for (unsigned int k = 0; k < toyparams_vector.at(j)->nparameters(); k++) {
        parameter* param = toyparams_vector.at(j)->get_parameter(k);
        RooRealVar* FinalParameter = (RooRealVar*)mc->GetParametersOfInterest()->find(param->get_name().c_str());
            double par_value = FinalParameter->getVal();
            FinalParameter->setVal(par_value-3*FinalParameter->getError());
            double nll1 = nll_hf->getVal();
            FinalParameter->setVal(par_value+3*FinalParameter->getError());
            double nll2 = nll_hf->getVal();
            FinalParameter->setVal(par_value);
            double approx_max = max(nll1,nll2);
            RooPlot *ParameterFrame = FinalParameter->frame("");
            nll_hf->plotOn(ParameterFrame, RooFit::ShiftToZero());
            ParameterFrame->SetMinimum(0.);
            ParameterFrame->SetMaximum(100.);//approx_max);
            ParameterFrame->Draw();
            gPad->SaveAs(TString("nll_")+param->get_name().c_str()+TString(".pdf"));
            ParameterFrame->Write();
        }
    }
    fnll->Close();
    
    
    output_simpleFit->Close();


    
    */

    char fitted_yields_prefix[256];
    strcpy(fitted_yields_prefix, "fitted_yields");
    saveResult(mc, common_parameters, toyparams_vector);
    saveFittedYields(cfg->outputFolder(), fitted_yields_prefix, name_suffix_vector, isocats, mc, common_parameters, toyparams_vector); 

    /*for (unsigned int i = 0; i<num_channels; i++){
      idx->setIndex(i);
      //if (cfg->ToyFit()==0) {
        plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,x_vector[i],cfg->outputFolder(), plotnames[0].c_str(), labels[0].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjectionInBins(idx, datahist_map_plotting, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjectionInBins(idx, datahist_map_plotting, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      }
    else {
        plotProjection(idx,toy_dataset_combined, (RooSimultaneous*) model,x_vector[i],cfg->outputFolder(), plotnames[0].c_str(), labels[0].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjection(idx,toy_dataset_combined, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjection(idx,toy_dataset_combined, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjectionInBins(idx,toy_dataset_combined, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
        plotProjectionInBins(idx,toy_dataset_combined, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      }
    }*/
  
    
    //model->Print("t");
  }

  else {

    LoadStartValues(mc, common_parameters, toyparams_vector);

    for (unsigned int i = 0; i<num_channels; i++){
      idx->setIndex(i);
      plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,x_vector[i],cfg->outputFolder(), plotnames[0].c_str(), labels[0].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      plotProjection(idx, datahist_map_plotting, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      plotProjectionInBins(idx, datahist_map_plotting, (RooSimultaneous*) model,y_vector[i],cfg->outputFolder(), plotnames[1].c_str(), labels[1].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
      plotProjectionInBins(idx, datahist_map_plotting, (RooSimultaneous*) model,z_vector[i],cfg->outputFolder(), plotnames[2].c_str(), labels[2].c_str(), mc,w, cfg->projections(), fin, name_suffix_vector[i].c_str(), error_histos[i]);
    }
  }
  
  return 0;
  
}
