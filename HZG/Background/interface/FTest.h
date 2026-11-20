#ifndef FLASHGGFINALFIT_HZG_F_TEST_H
#define FLASHGGFINALFIT_HZG_F_TEST_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdfCombine.h"
#include "RooRealVar.h"
#include "TString.h"

#include "BackgroundFitOptions.h"
#include "PdfModelBuilder.h"

struct FTestCategoryInputs {
  int cat = 0;
  int catOffset = 0;
  std::string catName;
  TString runPeriod;
  std::string energyExt;
  std::vector<std::string> functionClasses;
  RooRealVar* mass = nullptr;
  RooDataSet* data = nullptr;
  RooDataSet* dataSB = nullptr;
  RooDataSet* dataFull = nullptr;
  PdfModelBuilder* pdfsModel = nullptr;
  FILE* dfile = nullptr;
  FILE* resFile = nullptr;
};

struct FTestCategoryResult {
  FTestCategoryResult();

  std::map<std::string, int> choices;
  std::map<std::string, std::vector<int> > choicesEnvelope;
  std::map<std::string, RooAbsPdf*> pdfs;
  std::map<std::string, RooAbsPdf*> allPdfs;
  RooArgList storedPdfs;
  int simpleBestFitPdfIndex = 0;
  double minimumNLLSoFar = 1.e10;
};

FTestCategoryResult runEnvelopeFTest(const FTestCategoryInputs& input, BackgroundFitOptions& options);
int getBestFitFunction(RooMultiPdf* bkg, RooDataSet* data, RooDataSet* dataSB, RooCategory* cat, const BackgroundFitOptions& options, bool silent=false);

#endif
