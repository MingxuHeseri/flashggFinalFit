#ifndef FLASHGGFINALFIT_HZG_PLOT_UTILS_H
#define FLASHGGFINALFIT_HZG_PLOT_UTILS_H

#include <map>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdfCombine.h"
#include "RooRealVar.h"
#include "TString.h"

#include "BackgroundFitOptions.h"

void configureBackgroundPlotStyle();
void plot(RooRealVar* mass, RooAbsPdf* pdf, RooDataSet* data, std::string name, std::vector<std::string> flashggCats, int status, double* prob, BackgroundFitOptions& options);
void plot(RooRealVar* mass, RooMultiPdf* pdfs, RooCategory* catIndex, RooDataSet* data, std::string name, std::vector<std::string> flashggCats, int cat, int bestFitPdf, const BackgroundFitOptions& options);
void plot(RooRealVar* mass, std::map<std::string, RooAbsPdf*> pdfs, RooDataSet* data, std::string name, std::vector<std::string> flashggCats, int cat, int bestFitPdf, const BackgroundFitOptions& options);
void plot_corr(RooFitResult* modelFit, TString outputPath);

#endif
