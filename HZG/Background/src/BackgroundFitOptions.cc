#include "../interface/BackgroundFitOptions.h"

void BackgroundFitOptions::updateBinningFromMassRange() {
  nBinsForMass = 4 * static_cast<int>(mggHigh - mggLow);
}

RooRealVar& backgroundIntLumi() {
  static RooRealVar intLumi("IntLumi", "hacked int lumi", 1000.);
  return intLumi;
}

RooRealVar& backgroundSqrtS() {
  static RooRealVar sqrtS("SqrtS", "SqrtS", 13.);
  return sqrtS;
}
