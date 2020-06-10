#ifndef MULTIGAUSS_H
#define MULTIGAUSS_H
#include "RooGlobalFunc.h"
#include <stdlib.h>
#include "TMatrixDSym.h"
#include "RooMultiVarGaussian.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "TStopwatch.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MetropolisHastings.h"
#include "RooStats/MarkovChain.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/ProposalFunction.h"
#include "RooStats/PdfProposal.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/LikelihoodInterval.h"

#include <vector>
#include <map>
#include <unordered_map>

#include <string>

using namespace std;
using namespace RooFit;
using namespace RooStats;
namespace multifit {
    class MultiGauss{
    public:
        MultiGauss();
        ~MultiGauss();
        unordered_map<std::string, std::vector<Double_t> > defineValues();
        map<Int_t, std::string> stupidNaming();
        TMatrixDSym covMat(map<Int_t, std::string> stupidMap,unordered_map <std::string, std::vector<Double_t> > values);
        map<std::string, std::vector<Double_t> >  function();
    };
}
#endif