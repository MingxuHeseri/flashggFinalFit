#ifndef FLASHGGFINALFIT_HZG_SPURIOUS_SIGNAL_TEST_H
#define FLASHGGFINALFIT_HZG_SPURIOUS_SIGNAL_TEST_H

#include "RooRealVar.h"
#include "TString.h"

#include "BackgroundFitOptions.h"

bool SpuriousSignalTest(RooRealVar* inputMass, TString cat, TString funcType, int order, int sig, TString runPeriod, const BackgroundFitOptions& options);

#endif
