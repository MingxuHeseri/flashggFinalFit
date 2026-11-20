#ifndef FLASHGGFINALFIT_HZG_BACKGROUND_FIT_OPTIONS_H
#define FLASHGGFINALFIT_HZG_BACKGROUND_FIT_OPTIONS_H

#include <cmath>
#include <string>
#include <vector>

#include "RooRealVar.h"
#include "TRandom3.h"

struct BackgroundFitOptions {
  std::string fileName;
  int nCats = 5;
  int singleCategory = -1;
  int catOffset = 0;
  int nsigIn = 0;
  std::string datfile = "dat/fTest.dat";
  std::string outDir = "plots/fTest";
  std::string outfilename;
  std::string functionList = "BernsteinStepxGauFast,PowerLawStepxGau,ExponentialStepxGau,LaurentStepxGau,ModGaus,ExpModGauss";
  bool useDiGaussian = false;
  bool is2011 = false;
  bool verbose = false;
  bool saveMultiPdf = false;
  int isFlashgg = 1;
  bool isData = false;
  std::string flashggCatsStr = "UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag";
  std::vector<std::string> flashggCats;
  std::string year = "2016";

  bool blind = true;
  bool blindFit = false;
  bool doChi2Fit = false;
  bool plotCorr = false;
  bool runFtestCheckWithToys = false;
  bool skipSSTest = false;
  bool objectiveOnly = false;
  bool skipCrossObjectiveDiagnostics = false;
  double mggLow = 95.;
  double mggHigh = 170.;
  int nBinsForMass = 4 * (170 - 95);
  double blindLow = 120.;
  double blindHigh = 130.;
  int maxTries = 1;
  double minimizerEps = 0.001;
  int minimizerStrategy = 1;
  double upperEnvThreshold = 0.05;
  double minimumGofProbability = 0.10;
  double minimumKsProbability = 0.20;
  double gofBinsPerGeV = 1.;
  std::string sidebandChi2Mode = "integral";
  std::string gofIntegralMode = "auto";
  std::string plotNormalizationMode = "auto";
  int minPdfOrder = 1;
  int maxPdfOrder = 6;
  std::string candidateOrdersStr;
  std::vector<int> candidateOrders;
  int minimumCovQual = 3;
  double maximumEdm = 0.001;
  bool requireFitStatusZero = false;
  bool gateFitQuality = true;

  std::string sstInputPath;
  std::string sstSignalPath;

  TRandom3 randomGen;

  // Boundary proximity is diagnostic; it is not a numerical-quality veto.
  bool passesNumericalFitQuality(int status, int covarianceQuality, double edm) const {
    return (!requireFitStatusZero || status == 0) &&
           covarianceQuality >= minimumCovQual &&
           std::isfinite(edm) && edm <= maximumEdm;
  }

  void updateBinningFromMassRange();
};

RooRealVar& backgroundIntLumi();
RooRealVar& backgroundSqrtS();

#endif
