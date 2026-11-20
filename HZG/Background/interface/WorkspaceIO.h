#ifndef FLASHGGFINALFIT_HZG_WORKSPACE_IO_H
#define FLASHGGFINALFIT_HZG_WORKSPACE_IO_H

#include <string>

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TFile.h"

#include "BackgroundFitOptions.h"
#include "FTest.h"

void transferMacros(TFile* inFile, TFile* outFile);
bool writeMultiPdfForCategory(RooWorkspace* outputws, FTestCategoryResult& result, RooRealVar* mass, RooDataSet* data, RooDataSet* dataSB, RooDataSet* dataFull, int cat, const std::string& energyExt, BackgroundFitOptions& options, bool verbose);

#endif
