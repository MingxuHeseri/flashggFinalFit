#ifndef PdfModelBuilder_h 
#define PdfModelBuilder_h

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooProduct.h"
#include "RooWorkspace.h"
#include "TFile.h"

using namespace std;
using namespace RooFit;

class PdfModelBuilder {
  
  public:
    void setObsVar(RooRealVar *var);
    void setParameterConfigCategory(const string& cat);

    RooAbsPdf *getBernsteinStepxGauFast(string prefix, int order, bool doBlind, bool di_gauss=false);
    RooAbsPdf *getBernsteinErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getEmpiricalTail(string prefix, int order, bool doBlind, const string& family);
    RooAbsPdf *getExponentialErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getPowerLawErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getLaurentErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getExpBernsteinErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getExpChebychevErfTurnOn(string prefix, int order, bool doBlind);
    RooAbsPdf *getBernsteinStepxGauFastUnconstrained(string prefix, int order, bool doBlind);
    RooAbsPdf *getBernsteinStepxGauAnalytic(string prefix, int order, bool doBlind);
    RooAbsPdf *getTurnOnExpChebychev(string prefix, int order, bool doBlind);
    RooAbsPdf *getHZGErfExp(string prefix, bool doBlind);
    RooAbsPdf *getHZGPowerExp(string prefix, bool doBlind);
    RooAbsPdf *getHZGErfBernstein(string prefix, int order, bool doBlind);
    RooAbsPdf *getFkp(string prefix, int order, bool doBlind);
    RooAbsPdf *getExponentialStepxGau(string prefix, int order, bool doBlind, bool di_gauss, bool const_cp1);
    RooAbsPdf *getPowerLawStepxGau(string prefix, int order, bool doBlind, bool di_gauss, bool const_cp1);
    RooAbsPdf *getLaurentStepxGau(string prefix, int order, bool doBlind, bool di_gauss, bool const_cp1);
    RooAbsPdf *getModGaus(string prefix, bool doBlind);
    RooAbsPdf* getBernstein(string prefix, int order, bool doBlind,
                            const string& configFamily="Bernstein");
    RooAbsPdf* getChebychev(string prefix, int order);
    RooAbsPdf* getPowerLaw(string prefix, int order, bool doBlind);
    RooAbsPdf* getExponential(string prefix, int order, bool doBlind);
    RooAbsPdf* getLaurentSeries(string prefix, int order);
    RooAbsPdf* getExpModGaussian(string prefix, bool doBlind);
    RooAbsPdf* getAsymGenGaussian(string prefix, bool doBlind);

  private:
    RooRealVar* obs_var = nullptr;
    std::string parameter_config_category;
    map<string,RooRealVar*> params;
    map<string,RooFormulaVar*> prods;
    map<string,RooAbsPdf*> utilities;
};

#endif
