#ifndef FLASHGGFINALFIT_HZG_GOODNESS_OF_FIT_H
#define FLASHGGFINALFIT_HZG_GOODNESS_OF_FIT_H

#include <string>

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

#include "BackgroundFitOptions.h"

double getGoodnessOfFit(RooRealVar* mass, RooAbsPdf* pdf, RooDataSet* data, std::string name, BackgroundFitOptions& options, bool sidebandOnly=false, double* chi2Value=nullptr, int* ndofValue=nullptr);
double getKSProb(RooRealVar* mass, RooAbsPdf* pdf, RooDataSet* data, std::string name, const BackgroundFitOptions& options);

#endif
