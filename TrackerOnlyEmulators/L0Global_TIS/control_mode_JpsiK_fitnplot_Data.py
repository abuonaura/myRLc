import ROOT as r
import preselection
r.gROOT.ProcessLine('.x RooIpatia2.cxx++')
def main(data_file):
##################
# Observables

    JpsiKmass = r.RooRealVar("Bplus_M","M(J/#psi K)",(5280-50),(5280+50),'MeV/c^{2}')
    JpsiKmass.setBins(100)

##################
# Data

    print 'Loading data ...'

    data_file_full = r.TFile(data_file)
    tree_full = data_file_full.Get("DecayTree")
    tree_full.SetBranchStatus("*",0)
    for br in preselection.br_list: tree_full.SetBranchStatus(br,1)


    data_file_filtered = r.TFile(data_file[:-5]+"_filtered.root","recreate")
    tree_filtered = tree_full.CopyTree(preselection.cut_string + " && (Bplus_L0MuonDecision_TOS || Bplus_L0Global_TIS) && Bplus_Hlt1TrackMuonDecision_TOS")
    data_ =r.TH1F("data_","data_",JpsiKmass.getBins(),JpsiKmass.getMin(),JpsiKmass.getMax())
    tree_filtered.Draw("Bplus_M>>data_")
    data = r.RooDataHist("data","data",r.RooArgList(JpsiKmass),data_)

    print 'Data loaded.'

##################
# Model

    Bplus_mass = r.RooRealVar("Bplus_mass","#mu_{B^{+}}",5279.32,5250.,5300.)
    Bplus_width = r.RooRealVar("Bplus_width","#sigma_{B^{+}}",20,1.,30.)
    Bplus_lambda = r.RooRealVar("Bplus_lambda","lambda",-2.76668e+00,-5,-1)
    Bplus_zeta = r.RooRealVar("zeta","zeta",1e-10)
    Bplus_beta = r.RooRealVar("beta","beta",1e-10)
    Bplus_a1 = r.RooRealVar("Bplus_a1","a1",1.65026e+00,0.,5.)
    Bplus_n1 = r.RooRealVar("Bplus_n1","n1",2.68539e+00,0.,10.)
    Bplus_a2 = r.RooRealVar("Bplus_a2","a2",2.29820e+00,0.,5.)
    Bplus_n2 = r.RooRealVar("Bplus_n2","n2",3.57608e+00,0.,20.)

    Bplus_lambda.setVal(-3.39379529291)
    Bplus_a1.setVal(2.42605958582)
    Bplus_n1.setVal(2.65971449342)
    Bplus_a2.setVal(2.66726455948)
    Bplus_n2.setVal(2.86003122367)

    for var in [Bplus_lambda,Bplus_a1,Bplus_n1,Bplus_a2,Bplus_n2]: var.setConstant()

    Bplus_pdf = r.RooIpatia2("Bplus_pdf","Bplus_pdf",JpsiKmass,Bplus_lambda,Bplus_zeta,Bplus_beta,Bplus_width,Bplus_mass,Bplus_a1,Bplus_n1,Bplus_a2,Bplus_n2)
    Bplus_yield = r.RooRealVar("Bplus_yield","N_{B^{+}}",2e+05,0.,1e+07)

    bkg_coef = r.RooRealVar("bkg_coef","c_{bkg}",0,-0.1,0.1)
#bkg_pdf = RooPolynomial("bkg_pdf","bkg_pdf",JpsiKmass,RooArgList(bkg_coef))
    bkg_pdf = r.RooExponential("bkg_pdf","bkg_pdf",JpsiKmass,bkg_coef)
    bkg_yield = r.RooRealVar("bkg_yield","N_{bkg}",8e+04,0.,1.e+06)

    model = r.RooAddPdf("model","model",r.RooArgList(Bplus_pdf,bkg_pdf),r.RooArgList(Bplus_yield,bkg_yield))

##################
# Fit

    result = model.fitTo(data,r.RooFit.Extended(),r.RooFit.Save(),r.RooFit.Minos(1))
    result.Print()

##################
# Plot

    fr = JpsiKmass.frame()
    data.plotOn(fr,r.RooFit.Name("data"))
    model.plotOn(fr,r.RooFit.Name("model"))
    model.plotOn(fr,r.RooFit.Components("Bplus_pdf"),r.RooFit.Name("Bplus_comp"),r.RooFit.LineColor(r.kRed))
    model.plotOn(fr,r.RooFit.Components("bkg_pdf"),r.RooFit.Name("bkg_comp"),r.RooFit.LineColor(r.kGreen+2),r.RooFit.LineStyle(2))

    r.gStyle.SetLegendBorderSize(1)
    leg = r.TLegend(0.6,0.7,0.95,0.95)
    leg.AddEntry(fr.findObject("data"),"2015 data","epl")
    leg.AddEntry(fr.findObject("model"),"Full model","l")
    leg.AddEntry(fr.findObject("Bplus_comp"),"B^{+}->J/#psi K^{+}","l")
    leg.AddEntry(fr.findObject("bkg_comp"),"Background","l")

    c = r.TCanvas("c","c")
    fr.Draw()
    leg.Draw()
    c.Print('JpsiK_mass_plot'+data_file[-25:-21]+'.root')
    c.Print('JpsiK_mass_plot'+data_file[-25:-21]+'.pdf')

    params = model.getParameters(r.RooArgSet(JpsiKmass))
    #params.printLatex()

##################
# Compute sWeights

    data_unbinned = r.RooDataSet("data_unbinned","data_unbinned",tree_filtered,r.RooArgSet(JpsiKmass))

    for var in [Bplus_mass,Bplus_width,bkg_coef]: var.setConstant(1)
    r.RooMsgService.instance().setSilentMode(True)
    sData = r.RooStats.SPlot("sData","sData",data_unbinned,model,r.RooArgList(Bplus_yield,bkg_yield))

    sw_var = r.MyStruct()
    br_sw_var = tree_filtered.Branch('sweight', r.AddressOf(sw_var,'afloat'), 'sweight/F')
    Bplus_ETA = r.MyStruct()
    br_Bplus_ETA = tree_filtered.Branch('Bplus_ETA', r.AddressOf(Bplus_ETA,'afloat'), 'Bplus_ETA/F')
    Jpsi_ETA = r.MyStruct()
    br_Jpsi_ETA = tree_filtered.Branch('Jpsi_ETA', r.AddressOf(Jpsi_ETA,'afloat'), 'Jpsi_ETA/F')
    Kplus_ETA = r.MyStruct()
    br_Kplus_ETA = tree_filtered.Branch('Kplus_ETA', r.AddressOf(Kplus_ETA,'afloat'), 'Kplus_ETA/F')

    for i in range(tree_filtered.GetEntries()):
        sw_var.afloat = data_unbinned.get(i).getRealValue("Bplus_yield_sw")
        tree_filtered.GetEntry(i)
        Bplus_ETA.afloat = preselection.pseudorapidity(getattr(tree_filtered,"Bplus_P"),getattr(tree_filtered,"Bplus_PZ"))
        Jpsi_ETA.afloat = preselection.pseudorapidity(getattr(tree_filtered,"Jpsi_P"),getattr(tree_filtered,"Jpsi_PZ"))
        Kplus_ETA.afloat = preselection.pseudorapidity(getattr(tree_filtered,"Kplus_P"),getattr(tree_filtered,"Kplus_PZ"))
        for br in [br_sw_var,br_Bplus_ETA,br_Jpsi_ETA,br_Kplus_ETA]: br.Fill()

    tree_filtered.Write()
    data_file_filtered.Close()
