#!/bin/env python
from PIDPerfScripts.StartScreen import *
from PIDPerfScripts.DataFuncs import *
from PIDPerfScripts.Definitions import *
from PIDPerfScripts.Binning import *
from PIDPerfScripts.DataFuncs import *
from PIDPerfScripts.PerfCalcFuncs import *
from PIDPerfScripts.RunDictFuncs import *

import ROOT
import sys
import os.path
import itertools
import argparse
import yaml

tSize = .06

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(1)

def optimize_binning(particleName, 
                     binningVariable,
                     StrippingVersion,
                     polarity,
                     dataset,
                     pidCut,
                     variable_min,
                     variable_max,
                     additionalCut,
                     numberOfInitialBins,
                     startWithIsopopulatedBins,
                     mergeBelow,
                     mergeAbove,
                     minimumBinWidth,
                     delta,
                     nSigma,
                     binningSchemeName,
                     outputFile):

    year = StrippingVersion

    variable_names = {      "P":"_P",
                            "PT":"_PT",
                            "ETA":"_Eta",
                            "nTracks":"nTracks",
                            "nSPDHits":"nSPDHits",
                            "Brunel_P":"Brunel_P",
                            "Brunel_PT":"Brunel_PT",
                            "Brunel_ETA":"Brunel_ETA",
                        }

    # The number of inital bins is multiplicated by a factor of 100 for reading the dataset.
    # From this large set of bins the inital binning scheme with isopopulated bins or bins with the same sizes is calculated.
    # After this step the merging of the bins with the number of inital bins numberOfInitialBins is started.
    numberOfReadingBins = 100 * numberOfInitialBins

    var = variable_names[binningVariable]

    name = ""
    title = "" 

    if not "Brunel" in var:
        vName = var.replace("_","")
    else:
        vName = var

    vName = vName.replace("Eta", "ETA")

    draw = particleName + var
    if (var == "nTracks") or (var == "nSPDHits"):
        draw = var

    #Binning
    print("numberOfReadingBins: ", numberOfReadingBins)
    print("variable_min: ", variable_min)
    print("variable_max: ", variable_max)

    #Create histogram from full PID sample
    VarList = ROOT.RooArgList(dataset.Get_Param(binningVariable))
    hFTot = ROOT.TH1D(str("hFTot"+name), str("hFTot"+name), numberOfReadingBins, variable_min, variable_max)
    dataset.fillHistogram(hFTot, VarList,additionalCut)

    #Create histogram from PID sample with cut applied
    hFPas = ROOT.TH1D(str("hFPas"+name), str("hFTot"+name), numberOfReadingBins, variable_min, variable_max)
    dataset.fillHistogram(hFPas, VarList,pidCut+"&&("+additionalCut+")")
    hEff = getEfficiencyAutoBin(hFTot,
                                hFPas,
                                numberOfInitialBins,
                                nSigma,
                                delta,
                                variable_min,
                                variable_max,
                                var,
                                minimumBinWidth,
                                particleName+"_"+var+"_"+title,
                                startWithIsopopulatedBins,
                                mergeBelow,
                                mergeAbove)

    directory = os.path.dirname(outputFile)
    if not os.path.exists(directory) and not directory == "":
        os.makedirs(directory)

    if os.path.exists(outputFile):
        pyFile = open(outputFile, "a")
    else:
        pyFile = open(outputFile, "w")
        pyFile.write("from PIDPerfScripts.Binning import *\n")
        pyFile.write("from PIDPerfScripts.Definitions import *\n")

    pyFile.write("AddBinScheme('{}','{}', '{}', {}, {})\n".format(particleName, vName, binningSchemeName, variable_min, variable_max))
    for i in range(2,hEff.GetNbinsX()+1):
        pyFile.write("AddBinBoundary('{}', '{}', '{}', {})\n".format(particleName, vName, binningSchemeName, hEff.GetBinLowEdge(i)))
    pyFile.write("\n")
    pyFile.close()

    histogramFile = ROOT.TFile(outputFile.replace(".py", ".root"), "RECREATE")
    histogramFile.ls()
    #Plot
    hPC = ROOT.TH1D("hPC", "", 100, variable_min, variable_max)
    print(draw)
    dataset.fillHistogram(hPC, VarList)
    scale(hPC)  
    hEff.SetMinimum(0)
    hEff.SetMaximum(1.1)    
    if (var == "_P"):
        hEff.SetXTitle("#it{P} [MeV/c^{2}]")
    if (var == "_PT"):
        hEff.SetXTitle("#it{P_T} [MeV/c^{2}]")
    if (var == "_Eta"):
        hEff.SetXTitle("#it{#eta}")
    if (var == "nTracks"):
        hEff.SetXTitle("N_{tracks}")
    if (var == "nSPDHits"):
        hEff.SetXTitle("N_{SPD hits}")
    if (var == "Brunel_P"):
        hEff.SetXTitle("#it{Brunel P} [MeV/c^{2}]")
    if (var == "Brunel_PT"):
        hEff.SetXTitle("#it{Brunel P_T} [MeV/c^{2}]")
    if (var == "Brunel_ETA"):
        hEff.SetXTitle("#it{Brunel #eta}")
    if (var == "Brunel_nTracks"):
        hEff.SetXTitle("Brunel N_{tracks}")
    hEff.SetStats(False)  
    canvas = ROOT.TCanvas("", "", 1600, 900)
    hEff.SetFillColor(ROOT.kBlue)
    hEff.SetLineColor(ROOT.kBlue)
    hEff.SetMarkerColor(ROOT.kBlue)
    hEff.SetMarkerStyle(10)
    hEff.SetMarkerSize(1.5)
    hEff.Draw("E1")
    hPC.Draw("hist same")
    hPC.GetYaxis().SetTitle("Arbitrary Units")
    print(draw)
    outputFile = outputFile.replace(".py", ".pdf")
    canvas.Print(outputFile)
    hEff.Write()
    hPC.Write()
    histogramFile.Close()

def getEfficiencyAutoBin(hTot,
                         hPas,
                         nBins,
                         nSigma,
                         delta,
                         minimum_value,
                         maximum_value,
                         var = "",
                         minimumBinWidth = 0,
                         filename="",
                         startWithIsopopulatedBins=True,
                         mergeBelow=-1e9,
                         mergeAbove=1e9,
                        ):

        print("getEfficiencyAutoBin ",nBins," ",nSigma," ",delta," ",var," ",minimumBinWidth)
        print("bins hTot:",hTot.GetNbinsX())
        print("nBins:",nBins)
        
        if startWithIsopopulatedBins:
                bins = getIsoBinning(hPas, nBins)
        else:
                mins = np.zeros(nBins+1)
                for i in range(0,nBins+1):
                        mins[i] = hPas.GetXaxis().GetXmin() + (hPas.GetXaxis().GetXmax() - hPas.GetXaxis().GetXmin() ) / nBins * i
                bins = ROOT.RooBinning(nBins, mins)

        
        bins_array = bins.array()
        bins_number = bins.numBoundaries()
        bins_array.SetSize(bins_number)
        bins_array = np.array(list(bins_array))

        if not maximum_value in bins_array:
            bins_array = np.append(bins_array,maximum_value)
            bins_number = bins_number + 1
            print(bins_array)
            bins = ROOT.RooBinning(bins_number-1, bins_array)

        hRTot = getRebin(hTot, bins_number-1, bins_array)
        hRPas = getRebin(hPas, bins_number-1, bins_array)

        checkHistogram(hRPas, hRTot)
        hEff = efficiency(hRPas, hRTot, "")
        print("hEff bins:",hEff.GetNbinsX())
        
        print("xmin, xmax: ", hRTot.GetXaxis().GetXmin(), hRTot.GetXaxis().GetXmax())
        print("xmin, xmax: ", hRPas.GetXaxis().GetXmin(), hRPas.GetXaxis().GetXmax())

        loop = True
        print("Start merging bins:")
        while loop: 
                tmp = ROOT.RooBinning(bins)
                #print(bins.numBoundaries() - 1)
                #print(tmp.numBoundaries() -1)
                #print(hEff.GetNbinsX())

                i = 1
                while i < hEff.GetNbinsX(): 
                        merge = False
                        x1  = hEff.GetBinLowEdge(i)
                        x2  = hEff.GetBinLowEdge(i+1)

                        if x2 > mergeAbove:
                            merge = True

                        if x2 < mergeBelow:
                            merge = True

                        if ROOT.TMath.Abs(x2 - x1) < minimumBinWidth:
                                merge = True

                        e1  = hEff.GetBinContent(i)
                        eE1 = hEff.GetBinError(i)
                        e2  = hEff.GetBinContent(i+1)
                        eE2 = hEff.GetBinError(i+1)

                        if e1 == e2:
                                merge = True
                        if (e1 == 1) and (e2 == 0):
                                merge = True
                        if (e1 == 0) and (e2 == 1):
                                merge = True
                        if (e1 == 1) or (e1 == 0) or (e2 == 1) or (e2 == 0):
                                merge = True
                        
                        if not merge:
                                if eE1 / e1 > 1:
                                        merge = True
                                if eE2 / e2 > 1:
                                        merge = True
                        d = 0
                        s = 0

                        if not merge:
                                d = ROOT.TMath.Abs(e1 - e2)
                                s = d
                                if e1 != 0:
                                        d = d / e1
                                        s = s / ROOT.TMath.Sqrt(ROOT.TMath.Power(eE1, 2) + ROOT.TMath.Power(eE2, 2))
                                #Merging based on delta and nSigma
                                if delta is not None and nSigma is not None:
                                    if (d < delta) and (s < nSigma):
                                         merge = True
                                #Merging based on delta only
                                elif delta is not None:
                                    if d < delta:
                                        merge = True
                                #Merging based on nSigma only
                                else:
                                    if s < nSigma:
                                        merge = True

                        if merge:
                                #print(" - removing ",tmp.nearestBoundary(hEff.GetBinLowEdge(i+1)))
                                tmp.removeBoundary(tmp.nearestBoundary(hEff.GetBinLowEdge(i+1)))
                                #print("Num boundaries:")
                                #print(tmp.numBoundaries()-1)
                                #print(bins.numBoundaries()-1)
                                i = i + 1

                        i = i + 1

                #print(bins.numBoundaries() - 1," -> ",tmp.numBoundaries() - 1)

                if bins.numBoundaries() == tmp.numBoundaries():
                        loop = False

                bins = ROOT.RooBinning(tmp)

                bins_array = bins.array()
                bins_number = bins.numBoundaries()
                bins_array.SetSize(bins_number)
                bins_array = np.array(list(bins_array))

                print(bins_array)

                hRTot = getRebin(hTot, bins.numBoundaries() - 1, bins_array)
                hRPas = getRebin(hPas, bins.numBoundaries() - 1, bins_array)

                hEff = efficiency(hRPas, hRTot, "")

        return hEff

def getIsoBinning(hHisto, nBins):
        fraction = 1. / nBins
        hIntegrated = hHisto.GetCumulative(1, "").Clone("hIntegrated")
        hIntegrated.Scale(1. / hIntegrated.GetBinContent(hIntegrated.GetNbinsX()))
        bins = 1
        mins = np.zeros(100000)
        for i in range(1, hIntegrated.GetNbinsX()+1):
                if (hIntegrated.GetBinContent(i) > fraction * bins):
                        bins = bins + 1
                        mins[bins-1] = hIntegrated.GetBinLowEdge(i)
        bins = bins + 1
        mins[0] = hIntegrated.GetXaxis().GetXmin()
        mins[bins-1] = hIntegrated.GetXaxis().GetXmax()
        if (bins-1 != nBins):
                print("Wrong rebinning",nBins,bins-1)
                if bins-1 > nBins:
                        for i in range(1, ROOT.TMath.Abs(nBins - (bins-1))):
                                mins[bins-1-i-1] = mins[bins-1-i]
                                mins[bins-1-i] = 0
                else:
                        return
        return ROOT.RooBinning(nBins, mins)

def getRebin(hHisto, bins, mins=None):
    hHHisto = hHisto.Clone("hRebin")
    hHHisto.Sumw2()
    if mins is None:
            if (bins == -1):
                bins = hHHisto.GetNbinsX()
            return hHHisto.Rebin(bins)
    else:
        return hHHisto.Rebin(bins, "", mins)

def getEfficiencyWithNPasNTot(nPas, nTot, name=""):
        pwd = ""
        if name != "":
                pwd = ROOT.gDirectory.CurrentDirectory().GetPath()
                ROOT.gDirectory.mkdir(name)
                ROOT.gDirectory.cd(name)

        eff  = 1
        effE = 0
        isNeg = False
        if nPas < 0:
                nPas = nPas * -1
                isNeg = True

        hPas = ROOT.TH1D("", "", 1, 0, 1)
        hPas.SetBinContent(1, nPas)
        hPas.SetBinError(1, ROOT.TMath.Sqrt(nPas))
        hTot = ROOT.TH1D("", "", 1, 0, 1)
        hTot.SetBinContent(1, nTot)
        hTot.SetBinError(1, ROOT.TMath.Sqrt(nTot))
        tEff = ROOT.TEfficiency(hPas, hTot)
        tGraph = tEff.CreateGraph()
        if name != "":
                hPas.Write("hPas")
                hTot.Write("hTot")
                tEff.Write("tEff")
                tGraph.Write("tGraph")

        effEU = tEff.GetEfficiencyErrorUp(1)
        effEL = tEff.GetEfficiencyErrorLow(1)
        eff  = tEff.GetEfficiency(1)
        effE = (effEU + effEL) / 2.
        if isNeg:
                eff = eff * -1

        if name != "":
                print(name)
        print("NPas = ",hPas.GetBinContent(1))
        print("NTot = ",hTot.GetBinContent(1))
        print("Eff  = ",eff," +- ",effE,"  (+",effEU," -",effEL,")")
        if effE != 0:
                print(" (",effE / eff * 100,"%)")
        if name != "":
                ROOT.gDirectory.cd(pwd)
        return

def efficiency(hPas, hTot, name=""):
        print("Calculate Efficiency: name:{}".format(name))

        if name != "":
                pwd = ROOT.gDirectory.CurrentDirectory().GetPath()
                ROOT.gDirectory.mkdir(name)
                ROOT.gDirectory.cd(name)

        hHPas = hPas.Clone("hPas")
        hHTot = hTot.Clone("hTot")

        eff  = 1
        effE = 0

        nPas = hHPas.GetEntries()
        iPas = hHPas.Integral()

        nTot = hHTot.GetEntries()
        iTot = hHTot.Integral()

        wPas = iPas
        statPas = np.zeros(10)
        hHPas.GetStats(statPas)
        if ROOT.TMath.Abs(statPas[0] - statPas[1]) > 1e-5:
                print("*** WARNING: Pas histogram filled with weights")
                hIPas = getInt(hHPas).Clone("hIPas")
                hHPas = hIPas.Clone("hHPas")
                iPas = hHPas.Integral()

        wTot = iTot
        statTot = np.zeros(10)
        hHTot.GetStats(statTot)
        if (ROOT.TMath.Abs(statTot[0] - statTot[1]) > 1e-5):
                print("*** WARNING: Tot histogram filled with weights")
                hITot = getInt(hHTot).Clone("hITot")
                hHTot = hITot.Clone("hHTot")
                iTot = hHTot.Integral()

        checkHistogram(hHPas, hHTot)

        if not CompareHistograms(hHPas, hHTot):
            tEff = 0
        else:
            hN = hHPas
            hD = hHTot
            for i in range(1, hN.GetNbinsX()+1):
                if hN.GetBinContent(i) > hD.GetBinContent(i):
                    print("*** Num != Den: ",i," (",hN.GetBinLowEdge(i),"-",hN.GetBinLowEdge(i+1),")  ",hN.GetBinContent(i)," != ",hD.GetBinContent(i))
                    hN.SetBinContent(i, hD.GetBinContent(i))
                    hN.SetBinError(i, hD.GetBinError(i))
            tEff = ROOT.TEfficiency(hN, hD)

        tGraph = tEff.CreateGraph()

        hEff = hHPas.Clone("hEff")
        hEff.Reset()
        for i in range(1,hEff.GetNbinsX()+1):
                if (hHPas.GetBinContent(i) != 0) and (hHTot.GetBinContent(i) != 0):
                        hEff.SetBinContent(i, tEff.GetEfficiency(hEff.GetBin(i)))
                        hEff.SetBinError(i, (tEff.GetEfficiencyErrorUp(hEff.GetBin(i)) + tEff.GetEfficiencyErrorLow(hEff.GetBin(i))) / 2.)

        hEff.SetMinimum(0)
        hEff.SetMaximum(1.1)

        if name != "":
                hHPas.Write("hPas")
                hHTot.Write("hTot")
                tEff.Write("tEff")
                tGraph.Write("tGraph")
                hEff.Write("hEff")

        print("NPas = ",iPas," (",nPas,", ",wPas,")")
        print("NTot = ",iTot," (",nTot,", ",wTot,")")

        getEfficiencyWithNPasNTot(iPas, iTot)

        hEff.SetMinimum(0.7)
        hEff.SetMaximum(1.0)

        if name != "":
                ROOT.gDirectory.cd(pwd)

        return hEff

def CompareHistograms(hHisto1, hHisto2):
    if hHisto1.GetNbinsX() != hHisto2.GetNbinsX() \
        or hHisto1.GetXaxis().GetXmin() != hHisto2.GetXaxis().GetXmin() \
        or hHisto1.GetXaxis().GetXmax() != hHisto2.GetXaxis().GetXmax():
        print("Different Histograms (CompareHistograms)")
        return 0
    else:
        return 1

def getInt(hHisto):
        mins = np.zeros(hHisto.GetNbinsX()+1)

        for i in range(1, hHisto.GetNbinsX()+1):
                mins[i-1] = hHisto.GetBinLowEdge(i)
        mins[hHisto.GetNbinsX()] = hHisto.GetBinLowEdge(hHisto.GetNbinsX() + 1)

        hInt = ROOT.TH1D("", "", hHisto.GetNbinsX(), mins)
        hInt.SetXTitle(hHisto.GetXaxis().GetTitle())

        for i in range(1, hInt.GetNbinsX()+1):
                if ((hHisto.GetBinContent(i) - int(hHisto.GetBinContent(i)) ) < 0.5):
                        hInt.SetBinContent(i, int(hHisto.GetBinContent(i)) )
                else:
                        hInt.SetBinContent(i, int(hHisto.GetBinContent(i)) + 1)
        hInt.SetEntries(hHisto.GetEntries())
        return hInt

def checkHistogram(hPas, hTot):
    for i in range(1,hPas.GetNbinsX()+1):
        if hPas.GetBinContent(i) < 0:
            print("hPas bin content < 0")
            print("*** WARNING: ",hPas.GetName()," = ",hPas.GetBinContent(i)," +- ",hPas.GetBinError(i)," @ ",i)
            hPas.SetBinContent(i, 0)
            hPas.SetBinError(i, 0)
        if hTot.GetBinContent(i) < 0:
            print("hTot bin content < 0")
            print("*** WARNING: ",hTot.GetName()," = ",hTot.GetBinContent(i)," +- ",hTot.GetBinError(i)," @ ",i)
            hTot.SetBinContent(i, 0)
            hTot.SetBinError(i, 0)
        if hPas.GetBinContent(i) > hTot.GetBinContent(i):
            print("hPas bin content > hTot bin content")
            print("*** WARNING: ",hPas.GetName()," = ",hPas.GetBinContent(i)," +- ",hPas.GetBinError(i)," @ ",i)
            print("*** WARNING: ",hTot.GetName()," = ",hTot.GetBinContent(i)," +- ",hTot.GetBinError(i)," @ ",i)
            hPas.SetBinContent(i, hTot.GetBinContent(i))
            hPas.SetBinError(i, hTot.GetBinError(i))

def getMax(hHisto):
    iMax = -1
    max = -1e9
    for i in range(1, hHisto.GetNbinsX()+1):
        if hHisto.GetBinContent(i) != 0:
            max = ROOT.TMath.Max(hHisto.GetBinContent(i) + hHisto.GetBinError(i), max)
    return max

def getMin(hHisto):
    iMin = -1
    min = 1e9
    for i in range(1, hHisto.GetNbinsX()+1):
        if hHisto.GetBinContent(i) != 0:    
            min = ROOT.TMath.Min(hHisto.GetBinContent(i) - hHisto.GetBinError(i), min)
    return min

def scale(hHisto):
    iBin = 0
    max = -1e9
    for i in range(1, hHisto.GetNbinsX()+1):
        if (hHisto.GetBinContent(i) > max):
            max = hHisto.GetBinContent(i)
            iBin = i
    norm = 1. / hHisto.GetBinContent(iBin)
    hHisto.SetMinimum(0)
    hHisto.SetMaximum(1.1)
    normE = 0
    entries = hHisto.GetEntries()
    for i in range(0,hHisto.GetNbinsX()+1):
        hHisto.SetBinError(i, ROOT.TMath.Sqrt(ROOT.TMath.Power(norm * hHisto.GetBinError(i), 2) + ROOT.TMath.Power(hHisto.GetBinContent(i) * normE, 2)))
        hHisto.SetBinContent(i, hHisto.GetBinContent(i) * norm)
    hHisto.SetEntries(entries)

if '__main__' == __name__:

    #read and parse configuration file
    argp = argparse.ArgumentParser()
    argp.add_argument( "configFile", metavar="INFILE", type=str,
        help="Provide input configuration file to parse")
    args = argp.parse_args()
    with open(str(args.configFile), 'r') as ymlfile:
        cfg = yaml.load(ymlfile)

    # set the sample version (stripping version in Run I, and Turbo version in Run II. Internally called StripVersion)
    StripVersion=cfg['sampleVersion']
    CheckStripVer(StripVersion)
    
    # set the magnet polarity
    MagPolarity=cfg['magnetPolarity']
    CheckMagPol(MagPolarity)

    # set the particle name
    PartName=cfg['particleName']
    CheckPartType(PartName)
    
    #Additional sample check for pA/Ap samples
    if 'pA' in StripVersion or 'Ap' in StripVersion:
        CheckStripVerPartNameMagPol(StripVersion,PartName,MagPolarity)

    # set the PID cuts
    PIDCut = cfg['pidCut']
    if not CheckCuts(PIDCut,StripVersion):
        raise ValueError("Invalid PID cut: \"%s\""%PIDCut)

    if (len(cfg['priorCut'])>0):
        if isinstance(cfg['priorCut'],str):
        	if not CheckCuts(cfg['priorCut'],StripVersion):
        		print("Invalid cut string %s" %str(cfg['priorCut']))
        elif isinstance(cfg['priorCut'],list):
            if not CheckCuts(cfg['priorCut'].join(" "),StripVersion):
            	print("Invalid cut string %s" %str(cfg['priorCut']))

 	# set correct names of cut variables
 	variables = DataSetVariables()
    if PartName == "e_B_Jpsi":
        particle = "e"
    else:
        particle = PartName

    # pid cuts
    rawCut = cfg['pidCut']
    for i in "()[]&!|><=:+-*/":
        rawCut = rawCut.replace(i," " + i +" ")
    cutted_list = [x for x in rawCut.split(" ") if x!=""]
    for i, cut in enumerate(cutted_list):
        if cut not in "()[]&!|><=:+-*/":
            try:
                float(cut)
            except ValueError:
                variable = variables[cut].getName(particle)
                cutted_list[i] = variable
    cutstringPID = ""
    for cut in cutted_list:
        cutstringPID += cut

    #additional cuts prior to pid cuts    
    rawCut = cfg['priorCut']
    for i in "()[]&!|><=:+-*/":
        rawCut = rawCut.replace(i," " + i +" ")
    cutted_list = [x for x in rawCut.split(" ") if x!=""]
    for i, cut in enumerate(cutted_list):
        if cut not in "()[]&!|><=:+-*/":
            try:
                float(cut)
            except ValueError:
                variable = variables[cut].getName(particle)
                cutted_list[i] = variable
    cutstringPrior= ""
    for cut in cutted_list:
        cutstringPrior += cut


    print(GetProtonPIDPartTypes())
    if PartName in GetProtonPIDPartTypes():
        cutstringPrior = cutstringPrior.replace(PartName, "P")
        cutstringPID = cutstringPID.replace(PartName, "P")
        
    #custom cuts probably need some more work in the core PIDCalib scripts, 
    #as additional branches are not in the standard variable list
    #customCuts = cfg['customCut']

    RunMin = cfg['minRun']
    RunMax = cfg['maxRun']
    maxFiles = cfg['maxFiles']

    if RunMin is not None:
        try:
            int(RunMin)
        except ValueError:
            print("Argument to --minRun ('%s') is not an integer'." %RunMin)
        if RunMax is None:
            print("Min run was specified as %s, but no max run was given." %RunMin)

    if RunMax is not None:
        try:
            int(RunMax)
        except ValueError:
            print("Argument to --maxRun ('%s') is not an integer'." %RunMax)
        if RunMin is None:
            print("Max run was specified as %s, but no min run was given." %RunMax)

    if maxFiles is not None:
        try:
            int(maxFiles)
        except ValueError:
            print("Argument to --maxFiles ('%s') is not an integer'." %maxFiles)
        if maxFiles is None:
            print("Max files was specified as %s, but no min run was given." %maxFiles)
   
    varName = cfg['varName']
    if varName=='':
        print("Argument to --varName is an empty string.")
    CheckVarName(varName)
    YVarName = ""
    ZVarName = ""

    if cfg['schemeName'] is None:
        cfg['schemeName'] = "mybinning"

    if cfg["outputFile"] is None:
        cfg["outputFile"] = "binning.py"

    if cfg['delta'] is None and cfg['nSigma'] is None:
        raise ValueError("Delta and nSigma are both None")
    #======================================================================
    # Print information
    #======================================================================
    print('========= Requested data samples ===========')
    print("Stripping version: %s" %StripVersion)
    print("Magnet polarity: %s" %MagPolarity)
    print("Particle name: %s" %PartName)
    if len(cfg['priorCut'])>0:
        print("Initial cuts: %s" %cfg['priorCut'])
    print("PID cut: %s" %cutstringPID)
    print("Corrected cuts: %s" %cutstringPrior)
    print(cfg)
    print('============================================')

    verbose = True
    minEntries = 0
    allowMissingDataSets = False

    #======================================================================
    # Loop over all calibration subsamples
    #======================================================================
    
    #Use GetFiles function for Run I data
    if "Turbo" not in StripVersion and "Electron" not in StripVersion:
        files = GetFiles(StripVersion,MagPolarity,PartName,RunMin,RunMax,maxFiles,verbose)
    elif "Turbo" in StripVersion:
        files = GetWGPFiles(StripVersion,MagPolarity,verbose)
    elif "Electron" in StripVersion:
        files = GetElectronFiles(StripVersion,MagPolarity,verbose)
    
    data_total = ROOT.RooDataSet()
    i = 0
    for file in files[:maxFiles]:
        DataSet_temp = GetDataSet(StripVersion, MagPolarity, PartName, cfg['priorCut'], cfg['pidCut'], varName, YVarName, ZVarName, file, verbose,
                           allowMissingDataSets, minEntries)
        if i==0:
            data_total = DataSet_temp
        else:
            data_total.append(DataSet_temp)
        i+=1

    optimize_binning(PartName, 
                        varName, 
                        StripVersion,
                        MagPolarity,
                        data_total,
                        cutstringPID,
                        cfg['minimum'],
                        cfg['maximum'],
                        cutstringPrior, 
                        cfg['numberOfInitialBins'], 
                        cfg['startWithIsopopulatedBins'],
                        cfg['mergeBelow'],
                        cfg['mergeAbove'],
                        cfg['minimumBinWidth'],
                        cfg['delta'],
                        cfg['nSigma'],
                        cfg["schemeName"],
                        cfg["outputFile"]
                        )
