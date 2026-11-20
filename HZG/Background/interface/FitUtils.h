#ifndef FLASHGGFINALFIT_HZG_FIT_UTILS_H
#define FLASHGGFINALFIT_HZG_FIT_UTILS_H

#include <string>

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "BackgroundFitOptions.h"

RooFitResult* Chi2Fit(RooAbsReal& chi2, int printLevel=-1, double eps=100, bool offset=false, int strategy=0);
void print_params(RooAbsPdf* pdf);
void runFit(RooAbsPdf* pdf, RooDataSet* data, RooDataSet* dataSB,
            double* nll, int* status, const BackgroundFitOptions& options,
            bool doBlind=false, const std::string& cat="",
            const std::string& typePrefix="", int order=-1, int nBins=-1,
            const RooAbsPdf* pdfToSave=nullptr, bool chi2fit=false,
            int* covQual=nullptr, double* edm=nullptr,
            int* boundaryHit=nullptr);

#endif
