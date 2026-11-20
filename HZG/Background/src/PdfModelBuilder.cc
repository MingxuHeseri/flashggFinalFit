#include "TCanvas.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "RooPlot.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooEffProd.h"
#include "RooNumIntConfig.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooErfProductPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooExpPolyErfPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooPowerErfPdf.h"
#include "RooProdPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooFitResult.h"
#include "RooRandom.h"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"

#include "../interface/PdfModelBuilder.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooDoubleCBFast.h"
#include "RooGaussModel.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"

#include "HiggsAnalysis/CombinedLimit/interface/HGGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZGThresholdPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooFkpPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

using namespace std;
using namespace RooFit;
using namespace boost;

namespace {
struct PdfParameterOverride {
  bool hasValue = false;
  bool hasMin = false;
  bool hasMax = false;
  bool hasConstant = false;
  double value = 0.;
  double min = 0.;
  double max = 0.;
  bool constant = false;
};

const boost::property_tree::ptree& pdfParameterConfig() {
  static boost::property_tree::ptree config;
  static bool loaded = false;
  if (loaded) return config;
  loaded = true;

  const char* envPath = std::getenv("HZG_BACKGROUND_PARAM_CONFIG");
  if (!envPath || std::string(envPath).empty()) {
    throw std::runtime_error(
        "HZG_BACKGROUND_PARAM_CONFIG must name the explicit PDF parameter JSON");
  }
  std::ifstream input(envPath);
  if (!input.good()) {
    throw std::runtime_error(
        std::string("Cannot read HZG_BACKGROUND_PARAM_CONFIG: ") + envPath);
  }
  try {
    boost::property_tree::read_json(envPath, config);
  } catch (const std::exception& error) {
    throw std::runtime_error(
        std::string("Invalid HZG_BACKGROUND_PARAM_CONFIG ") + envPath
        + ": " + error.what());
  }
  std::cout << "[INFO] Loaded PDF parameter configuration from "
            << envPath << std::endl;
  return config;
}

void mergeParameterOverride(const boost::property_tree::ptree& node,
                            PdfParameterOverride& target) {
  auto value = node.get_optional<double>("value");
  if (value) {
    target.hasValue = true;
    target.value = *value;
  }
  auto min = node.get_optional<double>("min");
  if (min) {
    target.hasMin = true;
    target.min = *min;
  }
  auto max = node.get_optional<double>("max");
  if (max) {
    target.hasMax = true;
    target.max = *max;
  }
  auto constant = node.get_optional<bool>("constant");
  if (constant) {
    target.hasConstant = true;
    target.constant = *constant;
  }
}

void applyPdfParameterConfig(const std::string& category,
                             const std::string& family,
                             int order,
                             const std::string& key,
                             RooRealVar* var) {
  if (!var) return;
  const boost::property_tree::ptree& config = pdfParameterConfig();
  if (config.empty()) return;

  PdfParameterOverride overrideConfig;
  bool found = false;

  auto mergePath = [&](const std::string& path) {
    if (auto node = config.get_child_optional(path)) {
      mergeParameterOverride(*node, overrideConfig);
      found = true;
    }
  };

  mergePath("global." + key);
  mergePath("pdfs." + family + ".params." + key);
  if (order >= 0) {
    mergePath("pdfs." + family + ".orders." + std::to_string(order) + "." + key);
  }
  mergePath("families." + family + ".params." + key);
  if (order >= 0) {
    mergePath("families." + family + ".orders." + std::to_string(order) + "." + key);
  }
  if (!category.empty()) {
    const std::string categoryPath = "categories." + category;
    mergePath(categoryPath + ".global." + key);
    mergePath(categoryPath + ".pdfs." + family + ".params." + key);
    if (order >= 0) {
      mergePath(categoryPath + ".pdfs." + family + ".orders." +
                std::to_string(order) + "." + key);
    }
    mergePath(categoryPath + ".families." + family + ".params." + key);
    if (order >= 0) {
      mergePath(categoryPath + ".families." + family + ".orders." + std::to_string(order) + "." + key);
    }
  }

  if (!found) return;
  if (overrideConfig.hasMin) var->setMin(overrideConfig.min);
  if (overrideConfig.hasMax) var->setMax(overrideConfig.max);
  if (overrideConfig.hasValue) var->setVal(overrideConfig.value);
  if (overrideConfig.hasConstant) var->setConstant(overrideConfig.constant);
}

// Keep each relative coefficient explicit. A small last component otherwise
// loses precision when reconstructed as one minus the other fractions.
RooArgList normalizedMixtureFractions(const std::string& prefix, int firstLabel,
                                      const RooArgList& coefficients,
                                      std::map<std::string, RooFormulaVar*>& functions) {
  RooArgList fractions;
  if (coefficients.getSize() == 1) return fractions;

  std::ostringstream expression;
  for (int i = 0; i < coefficients.getSize(); ++i) {
    if (i > 0) expression << "+";
    expression << "@" << i;
  }
  const std::string sumName = prefix + "Sum";
  auto* sum = new RooFormulaVar(sumName.c_str(), sumName.c_str(),
                                expression.str().c_str(), coefficients);
  functions.emplace(sumName, sum);
  for (int i = 0; i < coefficients.getSize(); ++i) {
    const std::string name = prefix + std::to_string(firstLabel + i);
    const RooArgList args(*coefficients.at(i), *sum);
    auto* fraction = new RooFormulaVar(name.c_str(), name.c_str(), "@0/@1", args);
    functions.emplace(name, fraction);
    fractions.add(*fraction);
  }
  return fractions;
}

RooAbsPdf* maskErfProduct(const std::string& prefix, RooRealVar& observable,
                          RooAbsPdf* model, bool doBlind,
                          const char* blindedTitle) {
  if (!doBlind) return model;

  RooFormulaVar* blind = new RooFormulaVar(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(observable));
  auto* masked = new RooEffProd(Form("%s_SB", prefix.c_str()), blindedTitle,
                                *model, *blind);
  if (dynamic_cast<RooPowerErfPdf*>(model)) {
    // Integrate the discontinuities of the 120--130 GeV mask adaptively.
    auto* config = masked->specialIntegratorConfig(true);
    config->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
    config->setEpsAbs(1e-14);
    config->setEpsRel(1e-10);
  }
  return masked;
}

RooAbsPdf* makeErfTurnOnProduct(const std::string& prefix,
                                RooRealVar& observable,
                                RooAbsReal& step,
                                RooAbsReal& sigma,
                                RooAbsPdf& tail,
                                const bool doBlind,
                                const char* blindedTitle) {
  RooAbsPdf* turnOn = new RooHZGErfTurnOnPdf(
      Form("%s_erfturnon", prefix.c_str()),
      Form("%s_erfturnon", prefix.c_str()), observable, step, sigma);
  // Both factors depend on the mass observable.  RooEffProd is RooFit's
  // dedicated PDF-times-efficiency implementation for this topology; using
  // RooProdPdf here builds an overlapping-observable RooProduct graph that
  // crashes ROOT 6.36's vectorized evaluator during fitTo().
  RooAbsPdf* model = new RooEffProd(
      prefix.c_str(), prefix.c_str(), tail, *turnOn);
  return maskErfProduct(prefix, observable, model, doBlind, blindedTitle);
}
}

void PdfModelBuilder::setObsVar(RooRealVar *var){
  obs_var=var;
}


void PdfModelBuilder::setParameterConfigCategory(const string& cat){
  parameter_config_category = cat;
}


RooAbsPdf* PdfModelBuilder::getChebychev(string prefix, int order){
  
  RooArgList *coeffList = new RooArgList();
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.01,-10.,10.);
    //RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    //prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*params[name]);
  }
  //RooChebychev *cheb = new RooChebychev(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  RooPolynomial *cheb = new RooPolynomial(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  return cheb;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(bern->GetName(),bern));

}


RooAbsPdf *PdfModelBuilder::getBernsteinStepxGauFast(string prefix, int order,
                                                      bool doBlind,
                                                      bool di_gauss){
  if (order < 2 || order >= 7) return nullptr;

  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()),
                                    Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()),
                                     Form("%s_sigma", prefix.c_str()), 4., 0.5, 15.);
  RooRealVar *sigma2 = new RooRealVar(Form("%s_sigma2", prefix.c_str()),
                                      Form("%s_sigma2", prefix.c_str()), 5., 0.1, 10.);
  RooRealVar *gauc = new RooRealVar(Form("%s_gauc", prefix.c_str()),
                                    Form("%s_gauc", prefix.c_str()), 0.79677, 0., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()),
                                    Form("%s_step", prefix.c_str()), 105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category, "BernsteinStepxGauFast", order, "mean", mean);
  applyPdfParameterConfig(parameter_config_category, "BernsteinStepxGauFast", order, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category, "BernsteinStepxGauFast", order, "sigma2", sigma2);
  applyPdfParameterConfig(parameter_config_category, "BernsteinStepxGauFast", order, "gauc", gauc);
  applyPdfParameterConfig(parameter_config_category, "BernsteinStepxGauFast", order, "step", step);

  RooArgList coeffList;
  for (int i = 0; i < order; ++i) {
    if (i == 0) {
      const string name = Form("%s_b0", prefix.c_str());
      RooRealVar *reference = new RooRealVar(name.c_str(), name.c_str(), 1.);
      reference->setConstant(true);
      params.insert(pair<string, RooRealVar*>(name, reference));
      coeffList.add(*reference);
      continue;
    }

    // A normalized Bernstein PDF is invariant under common coefficient
    // rescaling, so b0=1 defines the reference.  Parameterize every remaining
    // positive ratio in log space: direct non-negative coefficients create
    // hard zero boundaries and non-smooth covariance estimates when a small
    // basis weight is preferred.
    const string logName = Form("%s_logb%d", prefix.c_str(), i);
    RooRealVar *logRatio = new RooRealVar(
        logName.c_str(), logName.c_str(), std::log(0.1 * (i + 1)), -20., 5.);
    applyPdfParameterConfig(parameter_config_category,
                            "BernsteinStepxGauFast", order,
                            Form("logb%d", i), logRatio);
    const string ratioName = Form("%s_b%d", prefix.c_str(), i);
    RooFormulaVar *ratio = new RooFormulaVar(
        ratioName.c_str(), ratioName.c_str(), "exp(@0)", RooArgList(*logRatio));
    params.insert(pair<string, RooRealVar*>(logName, logRatio));
    prods.insert(pair<string, RooFormulaVar*>(ratioName, ratio));
    coeffList.add(*ratio);
  }

  // Keep the factorized step kernel inside the physical observable range.
  // Extending an unconstrained Bernstein polynomial into the FFT buffer can
  // create negative sampled values when the fitted endpoint coefficient is
  // negative, even though the old analytic convolution remains positive in
  // the fit range.
  RooAbsPdf *stepPdf = new RooGaussStepBernsteinFast(
      Form("%s_step_bernfast%d", prefix.c_str(), order),
      Form("%s_step_bernfast%d", prefix.c_str(), order),
      *obs_var, *step, coeffList, obs_var->getMax());

  RooGaussModel *gau = new RooGaussModel(
      Form("%s_gau_bernfast%d", prefix.c_str(), order),
      Form("%s_gau_bernfast%d", prefix.c_str(), order),
      *obs_var, *mean, *sigma);
  RooGaussModel *gau2 = new RooGaussModel(
      Form("%s_gau2_bernfast%d", prefix.c_str(), order),
      Form("%s_gau2_bernfast%d", prefix.c_str(), order),
      *obs_var, *mean, *sigma2);
  RooAbsPdf *resolution = gau;
  if (di_gauss) {
    gauc->setRange(0., 1.);
    resolution = new RooAddPdf(
        Form("%s_addgau_bernfast%d", prefix.c_str(), order),
        Form("%s_addgau_bernfast%d", prefix.c_str(), order),
        *gau2, *gau, *gauc);
  }

  RooFFTConvPdf *convolution = new RooFFTConvPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *resolution);
  convolution->setBufferFraction(0.5);

  if (doBlind) {
    RooGenericPdf *blindPdf = new RooGenericPdf(
        Form("%s_blind", prefix.c_str()),
        Form("%s_blind", prefix.c_str()),
        "((@0 < 120)? 1:((@0 > 130)? 1:0))",
        RooArgList(*obs_var));
    return new RooProdPdf(
        Form("%s_SB", prefix.c_str()),
        "Blinded PDF (exclude 120-130)",
        RooArgList(*convolution, *blindPdf));
  }
  return convolution;
}


RooAbsPdf *PdfModelBuilder::getBernsteinErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 2 || order > 6) return nullptr;

  // Keep the established HZG turn-on parameter names and ranges. In this
  // direct form, step is x0 and sigma is the Gaussian standard deviation.
  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 15.);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  // RooBernsteinFast fixes the first coefficient to one internally. Retain
  // the log-ratio convention of BernsteinStepxGauFast for every other
  // coefficient, so the smooth tail stays positive and identifiable.
  RooArgList coefficientRatios;
  for (int i = 1; i < order; ++i) {
    const string logName = Form("%s_logb%d", prefix.c_str(), i);
    RooRealVar *logRatio = new RooRealVar(
        logName.c_str(), logName.c_str(), std::log(0.1 * (i + 1)), -20., 5.);
    applyPdfParameterConfig(parameter_config_category,
                            "BernsteinErfTurnOn", order,
                            Form("logb%d", i), logRatio);
    const string ratioName = Form("%s_b%d", prefix.c_str(), i);
    RooFormulaVar *ratio = new RooFormulaVar(
        ratioName.c_str(), ratioName.c_str(), "exp(@0)", RooArgList(*logRatio));
    params.insert(pair<string, RooRealVar*>(logName, logRatio));
    prods.insert(pair<string, RooFormulaVar*>(ratioName, ratio));
    coefficientRatios.add(*ratio);
  }

  RooAbsPdf* model = new RooErfProductPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *step, *sigma,
      coefficientRatios, RooArgList(), 0);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                        "Blinded Bernstein erf-turn-on PDF");
}

// Positive, identifiable empirical tails; no turn-on or convolution nodes.


RooAbsPdf *PdfModelBuilder::getEmpiricalTail(
    string prefix, int order, bool doBlind, const string& family) {
  if (order < 1 || order > 3) return nullptr;
  if (family != "ExpPolyTail" && family != "PowerLogTail" &&
      family != "LaurentTail") return nullptr;
  RooArgList args(*obs_var);
  std::ostringstream formula;
  const string z = family == "PowerLogTail" ?
      "(log(@0/125.)/0.25)" : "((@0-125.)/35.)";
  if (family == "LaurentTail") formula << "pow(@0/125.,-4.)";
  else formula << "exp(";
  for (int i = 1; i <= order; ++i) {
    const string name = Form("%s_p%d", prefix.c_str(), i);
    const bool laurent = family == "LaurentTail";
    RooRealVar* p = new RooRealVar(name.c_str(), name.c_str(),
        laurent ? -1. : (i == 1 ? -1.5 : 0.), -10., 10.);
    applyPdfParameterConfig(parameter_config_category, family, order,
                            Form("p%d", i), p);
    params.insert(pair<string, RooRealVar*>(name, p));
    args.add(*p);
    if (laurent) {
      // A fixed unit reference weight removes the common scale degeneracy.
      const int power = i == 1 ? -5 : (i == 2 ? -3 : -6);
      formula << "+exp(@" << i << ")*pow(@0/125.," << power << ")";
    } else {
      if (i > 1) formula << "+";
      formula << "@" << i << "*pow(" << z << "," << i << ")";
    }
  }
  if (family != "LaurentTail") formula << ")";
  RooAbsPdf* tail = new RooGenericPdf(prefix.c_str(), prefix.c_str(),
                                     formula.str().c_str(), args);
  if (!doBlind) return tail;
  RooFormulaVar* blind = new RooFormulaVar(
      Form("%s_blind", prefix.c_str()), "sideband mask",
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooEffProd(Form("%s_SB", prefix.c_str()),
                        "Blinded empirical tail", *tail, *blind);
}


RooAbsPdf *PdfModelBuilder::getExponentialErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 1 || order > 5 || order % 2 == 0) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 10.);
  applyPdfParameterConfig(parameter_config_category,
                          "ExponentialErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "ExponentialErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList slopes;
  RooArgList coefficients;
  for (int index = 0; index < order; index += 2) {
    const int label = index + 1;
    const string slopeName = Form("%s_p%d", prefix.c_str(), label);
    RooRealVar *slope = new RooRealVar(
        slopeName.c_str(), slopeName.c_str(), -0.05 * index, -0.5, 0.);
    applyPdfParameterConfig(parameter_config_category,
                            "ExponentialErfTurnOn", order,
                            Form("p%d", label), slope);
    params.insert(pair<string, RooRealVar*>(slopeName, slope));
    slopes.add(*slope);

    const string coefficientName = Form("%s_cp%d", prefix.c_str(), label);
    RooRealVar *coefficient = new RooRealVar(
        coefficientName.c_str(), coefficientName.c_str(),
        index == 0 ? 1. : 0.9, 0., 1.);
    applyPdfParameterConfig(parameter_config_category,
                            "ExponentialErfTurnOn", order,
                            Form("cp%d", label), coefficient);
    if (index == 0) coefficient->setConstant(true);
    params.insert(pair<string, RooRealVar*>(coefficientName, coefficient));
    coefficients.add(*coefficient);
  }

  RooArgList fractions = normalizedMixtureFractions(prefix + "_expfrac", 1, coefficients, prods);
  // Fractions retain the original individually-normalized tail semantics.
  RooAbsPdf* model = new RooErfProductPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *step, *sigma,
      slopes, fractions, 1);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                        "Blinded exponential erf-turn-on PDF");
}


RooAbsPdf *PdfModelBuilder::getPowerLawErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 1 || order > 5 || order % 2 == 0) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 10.);
  applyPdfParameterConfig(parameter_config_category,
                          "PowerLawErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "PowerLawErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList powerParameters;
  RooArgList coefficients;
  for (int index = 0; index < order; index += 2) {
    const int label = index + 1;
    const string powerName = Form("%s_p%d", prefix.c_str(), label);
    RooRealVar *power = new RooRealVar(
        powerName.c_str(), powerName.c_str(), -4. - index, -15., 0.);
    applyPdfParameterConfig(parameter_config_category,
                            "PowerLawErfTurnOn", order,
                            Form("p%d", label), power);
    params.insert(pair<string, RooRealVar*>(powerName, power));
    powerParameters.add(*power);

    const string coefficientName = Form("%s_cp%d", prefix.c_str(), label);
    RooRealVar *coefficient = new RooRealVar(
        coefficientName.c_str(), coefficientName.c_str(),
        index == 0 ? 1. : 0.1 * label, 0., 1.);
    applyPdfParameterConfig(parameter_config_category,
                            "PowerLawErfTurnOn", order,
                            Form("cp%d", label), coefficient);
    if (index == 0) coefficient->setConstant(true);
    params.insert(pair<string, RooRealVar*>(coefficientName, coefficient));
    coefficients.add(*coefficient);
  }

  RooArgList fractions = normalizedMixtureFractions(prefix + "_powfrac", 1, coefficients, prods);
  RooAbsPdf* model = new RooPowerErfPdf(prefix.c_str(), prefix.c_str(),
      *obs_var, *step, *sigma, powerParameters, fractions);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                              "Blinded power-law erf-turn-on PDF");
}


RooAbsPdf *PdfModelBuilder::getLaurentErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 2 || order > 6) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      5., 0.5, 10.);
  applyPdfParameterConfig(parameter_config_category,
                          "LaurentErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "LaurentErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList powerParameters;
  RooArgList coefficients;
  for (int index = 0; index < order; ++index) {
    const string powerName = Form("%s_p%d", prefix.c_str(), index);
    RooRealVar *power = new RooRealVar(
        powerName.c_str(), powerName.c_str(), -4. - index);
    applyPdfParameterConfig(parameter_config_category,
                            "LaurentErfTurnOn", order,
                            Form("p%d", index), power);
    params.insert(pair<string, RooRealVar*>(powerName, power));
    powerParameters.add(*power);

    const string coefficientName = Form("%s_cp%d", prefix.c_str(), index);
    RooRealVar *coefficient = new RooRealVar(
        coefficientName.c_str(), coefficientName.c_str(),
        index == 0 ? 1. : std::max(0.1, 0.5 - 0.1 * index), 0., 1.);
    applyPdfParameterConfig(parameter_config_category,
                            "LaurentErfTurnOn", order,
                            Form("cp%d", index), coefficient);
    if (index == 0) coefficient->setConstant(true);
    params.insert(pair<string, RooRealVar*>(coefficientName, coefficient));
    coefficients.add(*coefficient);
  }

  RooArgList fractions = normalizedMixtureFractions(prefix + "_laufrac", 0, coefficients, prods);
  RooAbsPdf* model = new RooPowerErfPdf(prefix.c_str(), prefix.c_str(),
      *obs_var, *step, *sigma, powerParameters, fractions);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                              "Blinded Laurent erf-turn-on PDF");
}


RooAbsPdf *PdfModelBuilder::getExpBernsteinErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 1 || order > 6) return nullptr;

  const double low = obs_var->getMin();
  const double high = obs_var->getMax();
  const double width = high - low;
  if (!(width > 0.)) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 15.);
  applyPdfParameterConfig(parameter_config_category,
                          "ExpBernsteinErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "ExpBernsteinErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList coefficients;
  for (int index = 1; index <= order; ++index) {
    const string name = Form("%s_c%d", prefix.c_str(), index);
    RooRealVar *coefficient = new RooRealVar(name.c_str(), name.c_str(), 0., -10., 10.);
    applyPdfParameterConfig(parameter_config_category, "ExpBernsteinErfTurnOn", order, Form("c%d", index), coefficient);
    params.insert(pair<string, RooRealVar*>(name, coefficient));
    coefficients.add(*coefficient);
  }
  RooAbsPdf* model = new RooExpPolyErfPdf(prefix.c_str(), prefix.c_str(),
      *obs_var, *step, *sigma, coefficients, low, high, 0);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                       "Blinded exponential-polynomial erf-turn-on PDF");
}


RooAbsPdf *PdfModelBuilder::getExpChebychevErfTurnOn(
    string prefix, int order, bool doBlind) {
  if (order < 1 || order > 8) return nullptr;

  const double low = obs_var->getMin();
  const double high = obs_var->getMax();
  const double width = high - low;
  if (!(width > 0.)) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      low, low - 15., low + 15.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 20.);
  applyPdfParameterConfig(parameter_config_category,
                          "ExpChebychevErfTurnOn", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "ExpChebychevErfTurnOn", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList coefficients;
  for (int index = 1; index <= order; ++index) {
    const string name = Form("%s_c%d", prefix.c_str(), index);
    RooRealVar *coefficient = new RooRealVar(name.c_str(), name.c_str(), 0., -5., 5.);
    applyPdfParameterConfig(parameter_config_category, "ExpChebychevErfTurnOn", order, Form("c%d", index), coefficient);
    params.insert(pair<string, RooRealVar*>(name, coefficient));
    coefficients.add(*coefficient);
  }
  RooAbsPdf* model = new RooExpPolyErfPdf(prefix.c_str(), prefix.c_str(),
      *obs_var, *step, *sigma, coefficients, low, high, 1);
  return maskErfProduct(prefix, *obs_var, model, doBlind,
                       "Blinded exponential-polynomial erf-turn-on PDF");
}


RooAbsPdf *PdfModelBuilder::getBernsteinStepxGauFastUnconstrained(
    string prefix, int order, bool doBlind) {
  if (order < 2 || order > 6) return nullptr;

  // This is a validation family for the analytic Bernstein-step model: retain
  // its signed coefficient parameterization and fixed Gaussian mean, while
  // replacing only the analytic convolution by RooFFTConvPdf.
  RooRealVar *mean = new RooRealVar(
      Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  mean->setConstant(true);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 15.);
  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinStepxGauFastUnconstrained", order,
                          "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinStepxGauFastUnconstrained", order,
                          "step", step);
  params.insert(pair<string, RooRealVar*>(mean->GetName(), mean));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));

  RooArgList coefficients;
  for (int i = 0; i < order; ++i) {
    const string name = Form("%s_b%d", prefix.c_str(), i);
    RooRealVar *coefficient = nullptr;
    if (i == 0) {
      coefficient = new RooRealVar(name.c_str(), name.c_str(), 1.);
      coefficient->setConstant(true);
    } else {
      coefficient = new RooRealVar(name.c_str(), name.c_str(),
                                   0.1 * (i + 1), -15., 15.);
      applyPdfParameterConfig(parameter_config_category,
                              "BernsteinStepxGauFastUnconstrained", order,
                              Form("b%d", i), coefficient);
    }
    params.insert(pair<string, RooRealVar*>(name, coefficient));
    coefficients.add(*coefficient);
  }

  RooAbsPdf *stepPdf = new RooGaussStepBernsteinFast(
      Form("%s_step_bernfastunconstrained%d", prefix.c_str(), order),
      Form("%s_step_bernfastunconstrained%d", prefix.c_str(), order),
      *obs_var, *step, coefficients, obs_var->getMax());
  RooGaussModel *resolution = new RooGaussModel(
      Form("%s_gau_bernfastunconstrained%d", prefix.c_str(), order),
      Form("%s_gau_bernfastunconstrained%d", prefix.c_str(), order),
      *obs_var, *mean, *sigma);
  RooFFTConvPdf *convolution = new RooFFTConvPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *resolution);
  convolution->setBufferFraction(0.5);

  if (!doBlind) return convolution;
  RooGenericPdf *blindPdf = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()), "Blinded unconstrained fast Bernstein-step PDF",
      RooArgList(*convolution, *blindPdf));
}


RooAbsPdf *PdfModelBuilder::getBernsteinStepxGauAnalytic(
    string prefix, int order, bool doBlind) {
  if (order < 2 || order > 6) return nullptr;

  RooRealVar *mean = new RooRealVar(
      Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  mean->setConstant(true);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 15.);
  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinStepxGauAnalytic", order, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category,
                          "BernsteinStepxGauAnalytic", order, "step", step);
  params.insert(pair<string, RooRealVar*>(mean->GetName(), mean));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));

  RooArgList coefficients;
  for (int i = 0; i < order; ++i) {
    if (i == 0) {
      const string name = Form("%s_b0", prefix.c_str());
      RooRealVar *reference = new RooRealVar(name.c_str(), name.c_str(), 1.);
      reference->setConstant(true);
      params.insert(pair<string, RooRealVar*>(name, reference));
      coefficients.add(*reference);
      continue;
    }
    const string ratioName = Form("%s_b%d", prefix.c_str(), i);
    RooRealVar *ratio = new RooRealVar(
        ratioName.c_str(), ratioName.c_str(), 0.1 * (i + 1), -15., 15.);
    applyPdfParameterConfig(parameter_config_category,
                            "BernsteinStepxGauAnalytic", order,
                            Form("b%d", i), ratio);
    params.insert(pair<string, RooRealVar*>(ratioName, ratio));
    coefficients.add(*ratio);
  }

  RooAbsPdf *model = new RooGaussStepBernstein(
      prefix.c_str(), prefix.c_str(), *obs_var, *mean, *sigma, *step,
      coefficients);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()), "Blinded analytic Bernstein-step PDF",
      RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getTurnOnExpChebychev(
    string prefix, int order, bool doBlind) {
  if (order < 1 || order > 8) return nullptr;

  const double low = obs_var->getMin();
  const double high = obs_var->getMax();
  const double width = high - low;
  if (!(width > 0.)) return nullptr;

  RooRealVar *step = new RooRealVar(
      Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()),
      low, low - 15., low + 15.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 20.);
  applyPdfParameterConfig(parameter_config_category,
                          "TurnOnExpChebychev", order, "step", step);
  applyPdfParameterConfig(parameter_config_category,
                          "TurnOnExpChebychev", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(step->GetName(), step));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList formulaArgs(*obs_var, *step, *sigma);
  const string scaledMass = Form(
      "((2*(@0-%.17g)/%.17g)-1)", low, width);
  std::ostringstream exponent;
  for (int degree = 1; degree <= order; ++degree) {
    const string name = Form("%s_c%d", prefix.c_str(), degree);
    RooRealVar *coefficient = new RooRealVar(
        name.c_str(), name.c_str(), 0., -5., 5.);
    applyPdfParameterConfig(parameter_config_category,
                            "TurnOnExpChebychev", order,
                            Form("c%d", degree), coefficient);
    params.insert(pair<string, RooRealVar*>(name, coefficient));
    formulaArgs.add(*coefficient);
    if (degree > 1) exponent << "+";
    exponent << "@" << (degree + 2) << "*cos(" << degree
             << "*acos(" << scaledMass << "))";
  }
  const string formula =
      "0.5*(1+erf((@0-@1)/@2))*exp(" + exponent.str() + ")";
  RooAbsPdf *model = new RooGenericPdf(
      prefix.c_str(), prefix.c_str(), formula.c_str(), formulaArgs);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()),
      "Blinded turn-on exponential Chebyshev PDF", RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getHZGErfExp(string prefix, bool doBlind) {
  const double low = obs_var->getMin();
  RooRealVar *x0 = new RooRealVar(
      Form("%s_x0", prefix.c_str()), Form("%s_x0", prefix.c_str()),
      low, low - 15., low + 15.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 20.);
  RooRealVar *lambda = new RooRealVar(
      Form("%s_lambda", prefix.c_str()), Form("%s_lambda", prefix.c_str()),
      0.03, 1.e-5, 0.5);
  applyPdfParameterConfig(parameter_config_category, "HZGErfExp", 1, "x0", x0);
  applyPdfParameterConfig(parameter_config_category, "HZGErfExp", 1, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category, "HZGErfExp", 1, "lambda", lambda);
  params.insert(pair<string, RooRealVar*>(x0->GetName(), x0));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));
  params.insert(pair<string, RooRealVar*>(lambda->GetName(), lambda));

  RooAbsPdf *model = new RooHZGErfExpPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *x0, *sigma, *lambda);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()), "Blinded HZG erf-exponential PDF",
      RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getHZGPowerExp(string prefix, bool doBlind) {
  const double low = obs_var->getMin();
  RooRealVar *x0 = new RooRealVar(
      Form("%s_x0", prefix.c_str()), Form("%s_x0", prefix.c_str()),
      low - 5., low - 30., low - 0.1);
  RooRealVar *alpha = new RooRealVar(
      Form("%s_alpha", prefix.c_str()), Form("%s_alpha", prefix.c_str()),
      1., 0.01, 20.);
  RooRealVar *lambda = new RooRealVar(
      Form("%s_lambda", prefix.c_str()), Form("%s_lambda", prefix.c_str()),
      0.03, 1.e-5, 0.5);
  applyPdfParameterConfig(parameter_config_category, "HZGPowerExp", 1, "x0", x0);
  applyPdfParameterConfig(parameter_config_category, "HZGPowerExp", 1, "alpha", alpha);
  applyPdfParameterConfig(parameter_config_category, "HZGPowerExp", 1, "lambda", lambda);
  params.insert(pair<string, RooRealVar*>(x0->GetName(), x0));
  params.insert(pair<string, RooRealVar*>(alpha->GetName(), alpha));
  params.insert(pair<string, RooRealVar*>(lambda->GetName(), lambda));

  RooAbsPdf *model = new RooHZGPowerExpPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *x0, *alpha, *lambda);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()), "Blinded HZG power-exponential PDF",
      RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getHZGErfBernstein(
    string prefix, int order, bool doBlind) {
  if (order != 2 && order != 3) return nullptr;

  const double low = obs_var->getMin();
  RooRealVar *x0 = new RooRealVar(
      Form("%s_x0", prefix.c_str()), Form("%s_x0", prefix.c_str()),
      low, low - 15., low + 15.);
  RooRealVar *sigma = new RooRealVar(
      Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()),
      4., 0.5, 20.);
  applyPdfParameterConfig(parameter_config_category,
                          "HZGErfBernstein", order, "x0", x0);
  applyPdfParameterConfig(parameter_config_category,
                          "HZGErfBernstein", order, "sigma", sigma);
  params.insert(pair<string, RooRealVar*>(x0->GetName(), x0));
  params.insert(pair<string, RooRealVar*>(sigma->GetName(), sigma));

  RooArgList coefficients;
  RooRealVar *reference = new RooRealVar(
      Form("%s_b0", prefix.c_str()), Form("%s_b0", prefix.c_str()), 1.);
  reference->setConstant(true);
  params.insert(pair<string, RooRealVar*>(reference->GetName(), reference));
  coefficients.add(*reference);
  for (int index = 1; index <= order; ++index) {
    const string logName = Form("%s_logb%d", prefix.c_str(), index);
    RooRealVar *logRatio = new RooRealVar(
        logName.c_str(), logName.c_str(), std::log(0.1 * (index + 1)),
        -20., 5.);
    applyPdfParameterConfig(parameter_config_category, "HZGErfBernstein",
                            order, Form("logb%d", index), logRatio);
    const string ratioName = Form("%s_b%d", prefix.c_str(), index);
    RooFormulaVar *ratio = new RooFormulaVar(
        ratioName.c_str(), ratioName.c_str(), "exp(@0)", RooArgList(*logRatio));
    params.insert(pair<string, RooRealVar*>(logName, logRatio));
    prods.insert(pair<string, RooFormulaVar*>(ratioName, ratio));
    coefficients.add(*ratio);
  }

  RooAbsPdf *model = new RooErfBernsteinPdf(
      prefix.c_str(), prefix.c_str(), *obs_var, *x0, *sigma, coefficients);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(
      Form("%s_SB", prefix.c_str()), "Blinded HZG erf-Bernstein PDF",
      RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getFkp(string prefix, int order, bool doBlind) {
  if (order < 0 || order > 2) return nullptr;

  const double high = obs_var->getMax();
  RooRealVar *scale = new RooRealVar(
      Form("%s_scale", prefix.c_str()), Form("%s_scale", prefix.c_str()),
      high + 10., high + 0.1, 1000.);
  RooRealVar *b = new RooRealVar(
      Form("%s_b", prefix.c_str()), Form("%s_b", prefix.c_str()),
      4., 0., 100.);
  RooRealVar *a0 = new RooRealVar(
      Form("%s_a0", prefix.c_str()), Form("%s_a0", prefix.c_str()),
      7.2, 0., 100.);
  applyPdfParameterConfig(parameter_config_category, "Fkp", order, "scale", scale);
  applyPdfParameterConfig(parameter_config_category, "Fkp", order, "b", b);
  applyPdfParameterConfig(parameter_config_category, "Fkp", order, "a0", a0);
  params.insert(pair<string, RooRealVar*>(scale->GetName(), scale));
  params.insert(pair<string, RooRealVar*>(b->GetName(), b));
  params.insert(pair<string, RooRealVar*>(a0->GetName(), a0));

  RooAbsReal *a1 = nullptr;
  if (order >= 1) {
    RooRealVar *parameter = new RooRealVar(
        Form("%s_a1", prefix.c_str()), Form("%s_a1", prefix.c_str()),
        0., -20., 20.);
    applyPdfParameterConfig(parameter_config_category, "Fkp", order, "a1",
                            parameter);
    params.insert(pair<string, RooRealVar*>(parameter->GetName(), parameter));
    a1 = parameter;
  } else {
    a1 = new RooConstVar(Form("%s_a1", prefix.c_str()),
                         Form("%s_a1", prefix.c_str()), 0.);
  }

  RooAbsReal *a2 = nullptr;
  if (order >= 2) {
    RooRealVar *parameter = new RooRealVar(
        Form("%s_a2", prefix.c_str()), Form("%s_a2", prefix.c_str()),
        0., -20., 20.);
    applyPdfParameterConfig(parameter_config_category, "Fkp", order, "a2",
                            parameter);
    params.insert(pair<string, RooRealVar*>(parameter->GetName(), parameter));
    a2 = parameter;
  } else {
    a2 = new RooConstVar(Form("%s_a2", prefix.c_str()),
                         Form("%s_a2", prefix.c_str()), 0.);
  }

  RooAbsPdf *model = new RooFkpPdf(prefix.c_str(), prefix.c_str(), order,
                                    *obs_var, *scale, *b, *a0, *a1, *a2);
  if (!doBlind) return model;
  RooGenericPdf *blind = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(Form("%s_SB", prefix.c_str()), "Blinded Fk PDF",
                        RooArgList(*model, *blind));
}


RooAbsPdf *PdfModelBuilder::getPowerLawStepxGau(string prefix, int order, bool doBlind, bool di_gauss = false, bool const_cp1 = false){
  if (order%2==0 || order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 4., 0.5, 10.);
  // sigma->setConstant(true);
  RooRealVar *sigma2 = new RooRealVar(Form("%s_sigma2", prefix.c_str()), Form("%s_sigma2", prefix.c_str()), 5., 0.1, 10.);
  RooRealVar *gauc = new RooRealVar(Form("%s_gauc", prefix.c_str()), Form("%s_gauc", prefix.c_str()), 0.5, 0., 10.);
  // RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 97., 95., 105.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, "mean", mean);
  applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, "sigma2", sigma2);
  applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, "gauc", gauc);
  applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, "step", step);
  RooArgList powerList;
  RooArgList coefList;

  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1);
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -4.0-i, -15.0, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.1 * (i + 1), 0., 1.);
    if (i == 0) {
      if (order == 1) cp->setVal(1.0);
      else cp->setVal(0.1);
      if ((const_cp1)) cp->setConstant(true);
      // if ((const_cp1) && (order != 1)) cp->setConstant(true);
    }
    applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, Form("p%d", i + 1), p);
    applyPdfParameterConfig(parameter_config_category, "PowerLawStepxGau", order, Form("cp%d", i + 1), cp);
    powerList.add(*p);
    coefList.add(*cp);
  }
  RooAbsPdf *stepPdf = new RooGaussStepPowerLaw(
      Form("%s_step_pow%d", prefix.c_str(), order),
      Form("%s_step_pow%d", prefix.c_str(), order),
      *obs_var, *step, powerList, coefList, obs_var->getMax());
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_pow%d", prefix.c_str(), order), Form("%s_gau_pow%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooGaussModel *gau2 = new RooGaussModel(Form("%s_gau2_pow%d", prefix.c_str(), order), Form("%s_gau2_pow%d", prefix.c_str(), order), *obs_var, *mean, *sigma2);
  RooAddPdf *addgau = new RooAddPdf(Form("%s_addgau_pow%d", prefix.c_str(), order), Form("%s_addgau_pow%d", prefix.c_str(), order), *gau2, *gau, *gauc);
  string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
  RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), formula_SB.c_str(), RooArgList(*obs_var));

  if (di_gauss){
    RooFFTConvPdf *digauxpow = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *addgau);
    digauxpow->setBufferFraction(0.5);
    if (doBlind) {      
      RooProdPdf* digauxpow_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*digauxpow, *blindPdf)
      );
      return digauxpow_SB;
    }
    return digauxpow;
  } else {
    RooFFTConvPdf *gauxpow = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *gau);
    gauxpow->setBufferFraction(0.5);
    if (doBlind) {
      RooProdPdf* gauxpow_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*gauxpow, *blindPdf)
      );
      return gauxpow_SB;
    }
    return gauxpow;
  }
}


RooAbsPdf *PdfModelBuilder::getExponentialStepxGau(string prefix, int order, bool doBlind, bool di_gauss = false, bool const_cp1 = false){
  if (order%2==0 || order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 4., 0.5, 10.);
  // sigma->setConstant(true);
  RooRealVar *sigma2 = new RooRealVar(Form("%s_sigma2", prefix.c_str()), Form("%s_sigma2", prefix.c_str()), 5., 0.1, 10.);
  RooRealVar *gauc = new RooRealVar(Form("%s_gauc", prefix.c_str()), Form("%s_gauc", prefix.c_str()), 0.79677, 0., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, "mean", mean);
  applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, "sigma2", sigma2);
  applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, "gauc", gauc);
  applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, "step", step);
  RooArgList slopeList;
  RooArgList coefList;

  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1);
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -0.05*i, -0.5, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.9, 0., 1.);
    if (i == 0) {
      if (order == 1) cp->setVal(1.0);
      else cp->setVal(0.1);
      if ((const_cp1) && (order != 1)) cp->setConstant(true);
    }
    applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, Form("p%d", i + 1), p);
    applyPdfParameterConfig(parameter_config_category, "ExponentialStepxGau", order, Form("cp%d", i + 1), cp);
    slopeList.add(*p);
    coefList.add(*cp);
  }
  RooAbsPdf *stepPdf = new RooGaussStepExponential(
      Form("%s_step_exp%d", prefix.c_str(), order),
      Form("%s_step_exp%d", prefix.c_str(), order),
      *obs_var, *step, slopeList, coefList, obs_var->getMax());
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_exp%d", prefix.c_str(), order), Form("%s_gau_exp%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooGaussModel *gau2 = new RooGaussModel(Form("%s_gau2_exp%d", prefix.c_str(), order), Form("%s_gau2_exp%d", prefix.c_str(), order), *obs_var, *mean, *sigma2);
  RooAddPdf *addgau = new RooAddPdf(Form("%s_addgau_exp%d", prefix.c_str(), order), Form("%s_addgau_exp%d", prefix.c_str(), order), *gau2, *gau, *gauc);
  string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
  RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), formula_SB.c_str(), RooArgList(*obs_var));
  
  if (di_gauss){
    RooFFTConvPdf *digauxexp = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *addgau);
    if (doBlind) {      
      RooProdPdf* digauxexp_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*digauxexp, *blindPdf)
      );
      return digauxexp_SB;
    }
    return digauxexp;
  } else {
    RooFFTConvPdf *gauxexp = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *gau);
    if (doBlind) {
      RooProdPdf* gauxexp_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*gauxexp, *blindPdf)
      );
      return gauxexp_SB;
    }
    return gauxexp;
  }
}


RooAbsPdf *PdfModelBuilder::getLaurentStepxGau(string prefix, int order, bool doBlind, bool di_gauss = false, bool const_cp1 = false){
  if (order > 7) return NULL;
  RooRealVar *mean = new RooRealVar(Form("%s_mean", prefix.c_str()), Form("%s_mean", prefix.c_str()), 0.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma", prefix.c_str()), Form("%s_sigma", prefix.c_str()), 5., 0.5, 10.);
  RooRealVar *sigma2 = new RooRealVar(Form("%s_sigma2", prefix.c_str()), Form("%s_sigma2", prefix.c_str()), 1., 0.1, 10.);
  RooRealVar *gauc = new RooRealVar(Form("%s_gauc", prefix.c_str()), Form("%s_gauc", prefix.c_str()), 0.79677, 0., 10.);
  RooRealVar *step = new RooRealVar(Form("%s_step", prefix.c_str()), Form("%s_step", prefix.c_str()), 105., 95., 115.);
  applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, "mean", mean);
  applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, "sigma", sigma);
  applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, "sigma2", sigma2);
  applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, "gauc", gauc);
  applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, "step", step);
  RooArgList powerList;
  RooArgList coefList;

  for (int i = 0; i < order; i++){
    RooRealVar *cp = new RooRealVar(Form("%s_cp%d", prefix.c_str(), i), Form("%s_cp%d", prefix.c_str(), i), 0.5-0.1*i, 0, 1.0);
    RooRealVar *p = new RooRealVar(Form("%s_p%d", prefix.c_str(), i), Form("%s_p%d", prefix.c_str(), i), -4-i);
    if (i == 0) {
      if (order == 1) cp->setVal(1.0);
      else cp->setVal(2.0);
      if ((const_cp1) && (order != 1)) cp->setConstant(true);
    }
    applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, Form("p%d", i), p);
    applyPdfParameterConfig(parameter_config_category, "LaurentStepxGau", order, Form("cp%d", i), cp);
    powerList.add(*p);
    coefList.add(*cp);
  }

  RooAbsPdf *stepPdf = new RooGaussStepLaurent(
      Form("%s_step_lau%d", prefix.c_str(), order),
      Form("%s_step_lau%d", prefix.c_str(), order),
      *obs_var, *step, powerList, coefList, obs_var->getMax());
  RooGaussModel *gau = new RooGaussModel(Form("%s_gau_lau%d", prefix.c_str(), order), Form("%s_gau_lau%d", prefix.c_str(), order), *obs_var, *mean, *sigma);
  RooGaussModel *gau2 = new RooGaussModel(Form("%s_gau2_lau%d", prefix.c_str(), order), Form("%s_gau2_lau%d", prefix.c_str(), order), *obs_var, *mean, *sigma2);
  RooAddPdf *addgau = new RooAddPdf(Form("%s_addgau_lau%d", prefix.c_str(), order), Form("%s_addgau_lau%d", prefix.c_str(), order), *gau2, *gau, *gauc);
  string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
  RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), formula_SB.c_str(), RooArgList(*obs_var));
  
  if (di_gauss){
    RooFFTConvPdf *digauxlau = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *addgau);
    if (doBlind) {      
      RooProdPdf* digauxlau_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*digauxlau, *blindPdf)
      );
      return digauxlau_SB;
    }
    return digauxlau;
  } else {
    RooFFTConvPdf *gauxlau = new RooFFTConvPdf(prefix.c_str(), prefix.c_str(), *obs_var, *stepPdf, *gau);
    if (doBlind) {
      RooProdPdf* gauxlau_SB = new RooProdPdf(
          Form("%s_SB", prefix.c_str()),
          "Blinded PDF (exclude 120-130)",
          RooArgList(*gauxlau, *blindPdf)
      );
      return gauxlau_SB;
    }
    return gauxlau;
  }
}


RooAbsPdf *PdfModelBuilder::getModGaus(string prefix, bool doBlind) {
  // 定义 modified Gaussian 的参数
  RooRealVar *m0 = new RooRealVar(Form("%s_m0", prefix.c_str()), Form("%s_m0", prefix.c_str()), 105., 100., 130.);
  RooRealVar *vl = new RooRealVar(Form("%s_vl", prefix.c_str()), Form("%s_vl", prefix.c_str()), 2., 0., 15.);
  RooRealVar *vr = new RooRealVar(Form("%s_vr", prefix.c_str()), Form("%s_vr", prefix.c_str()), 1., -5., 5.);
  RooRealVar *s0 = new RooRealVar(Form("%s_s0", prefix.c_str()), Form("%s_s0", prefix.c_str()), 0.);
  RooRealVar *sl = new RooRealVar(Form("%s_sl", prefix.c_str()), Form("%s_sl", prefix.c_str()), 1., 1., 40.);
  RooRealVar *sh = new RooRealVar(Form("%s_sh", prefix.c_str()), Form("%s_sh", prefix.c_str()), 40., 1., 60.);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "m0", m0);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "vl", vl);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "vr", vr);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "s0", s0);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "sl", sl);
  applyPdfParameterConfig(parameter_config_category, "ModGaus", 1, "sh", sh);
  
  // 获取观测变量的范围用于归一化
  double xlow = obs_var->getMin();
  double xhigh = obs_var->getMax();
  
  // 创建 modified Gaussian PDF
  RooModGaus *modgaus = new RooModGaus(prefix.c_str(), prefix.c_str(),
                                       *obs_var, *m0, *vl, *vr, *s0, *sl,
                                       *sh, xlow, xhigh);
  
  // 盲化区域：排除 120-130 GeV（根据你的物理分析调整）
  string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
  RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), 
                                              formula_SB.c_str(), RooArgList(*obs_var));

  if (doBlind) {
    // 应用盲化：将 modified Gaussian 与盲化函数相乘
    RooProdPdf* modgaus_SB = new RooProdPdf(
        Form("%s_SB", prefix.c_str()),
        "Blinded Modified Gaussian (exclude 120-130)",
        RooArgList(*modgaus, *blindPdf)
    );
    return modgaus_SB;
  }
  
  return modgaus;
}


RooAbsPdf* PdfModelBuilder::getBernstein(string prefix, int order, bool doBlind,
                                        const string& configFamily){
  
  RooArgList *coeffList = new RooArgList();
  //coeffList->add(RooConst(1.0)); // no need for cnstant in this interface
  for (int i=0; i<order; i++){
    string name = Form("%s_p%d",prefix.c_str(),i);
    //params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),1.0,0.,5.)));
    RooRealVar *param = new RooRealVar(name.c_str(),name.c_str(),0.1*(i+1),-15.,15.);
    applyPdfParameterConfig(parameter_config_category, configFamily, order,
                            Form("p%d", i), param);
    RooFormulaVar *form = new RooFormulaVar(Form("%s_sq",name.c_str()),Form("%s_sq",name.c_str()),"@0*@0",RooArgList(*param));
    params.insert(pair<string,RooRealVar*>(name,param));
    prods.insert(pair<string,RooFormulaVar*>(name,form));
    coeffList->add(*prods[name]);
  }
  RooAbsPdf *bern = nullptr;
  if (order==1) {
	bern = new RooBernsteinFast<1>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else if (order==2) {
	bern = new RooBernsteinFast<2>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else if (order==3) {
	bern = new RooBernsteinFast<3>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else if (order==4) {
	bern = new RooBernsteinFast<4>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else if (order==5) {
	bern = new RooBernsteinFast<5>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else if (order==6) {
	bern = new RooBernsteinFast<6>(prefix.c_str(),prefix.c_str(),*obs_var,*coeffList);
  } else {
	return nullptr;
  }

  if (!doBlind) return bern;
  RooGenericPdf *blindPdf = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()),
      "((@0 < 120)? 1:((@0 > 130)? 1:0))", RooArgList(*obs_var));
  return new RooProdPdf(Form("%s_SB", prefix.c_str()),
                        "Blinded positive Bernstein PDF", RooArgList(*bern, *blindPdf));
}


RooAbsPdf *PdfModelBuilder::getPowerLaw(string prefix, int order, bool doBlind) {
  if (order % 2 == 0 || order > 7) return NULL;
  
  RooArgList formulaArgs;
  string formula = "(";
  
  formulaArgs.add(*obs_var);
  int xmin = obs_var->getMin();
  int xmax = obs_var->getMax();
  
  // 创建参数列表并构建公式
  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1);
    
    // 注意：这里去掉了原始代码中多余的order参数
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -4.0 - i, -15.0, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.1 * (i + 1), 0., 1.);
    
    // 设置第一个cp参数的初始值
    if (i == 0) {
      if (order == 1) {
        cp->setVal(1.0);
      } else {
        cp->setVal(0.1);
      }
    }
    
    formulaArgs.add(*p);
    formulaArgs.add(*cp);
    
    // 构建归一化的幂律函数公式
    // 原始公式: cp * x^p * (1+p) / (xmax^(1+p) - step^(1+p))
    // 修改后: cp * x^p / ∫ x^p dx (从xmin到xmax)
    if (i > 0) {
      formula += "+";
    }
    
    // 归一化的幂律函数公式
    // 对于p ≠ -1: ∫ x^p dx = (x^(p+1))/(p+1)
    // 归一化因子: 1/((xmax^(p+1) - xmin^(p+1))/(p+1))
    // 所以: cp * x^p / ((xmax^(p+1) - xmin^(p+1))/(p+1))
    formula += Form("@%d*TMath::Power(@0, @%d)/((TMath::Power(%d, 1+@%d)-TMath::Power(%d, 1+@%d))/(1+@%d))",
      formulaArgs.getSize() - 1,  // cp
      formulaArgs.getSize() - 2,  // p
      xmax,
      formulaArgs.getSize() - 2,  // p
      xmin,
      formulaArgs.getSize() - 2,  // p
      formulaArgs.getSize() - 2   // p
    );
  }
  
  formula += ")";
  cout << "PowerLaw PDF formula: " << formula << endl;
  
  // 创建幂律PDF
  RooGenericPdf *powPdf = new RooGenericPdf(
    Form("%s_pow%d", prefix.c_str(), order),
    Form("%s_pow%d", prefix.c_str(), order),
    formula.c_str(),
    formulaArgs
  );
  
  // 如果需要盲化处理
  if (doBlind) {
    string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
    RooGenericPdf *blindPdf = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()),
      Form("%s_blind", prefix.c_str()),
      formula_SB.c_str(),
      RooArgList(*obs_var)
    );
    
    RooProdPdf *powPdf_SB = new RooProdPdf(
      Form("%s_SB", prefix.c_str()),
      "Blinded PDF (exclude 120-130)",
      RooArgList(*powPdf, *blindPdf)
    );
    
    return powPdf_SB;
  }
  
  return powPdf;
}


RooAbsPdf *PdfModelBuilder::getExponential(string prefix, int order, bool doBlind) {
  if (order % 2 == 0 || order > 7) return NULL;
  
  RooArgList formulaArgs;
  string formula = "(";
  
  formulaArgs.add(*obs_var);
  int xmin = obs_var->getMin();
  int xmax = obs_var->getMax();
  
  // 创建参数列表并构建公式
  for (int i = 0; i < order; i += 2) {
    string pName = Form("%s_p%d", prefix.c_str(), i + 1);
    string cpName = Form("%s_cp%d", prefix.c_str(), i + 1);
    
    RooRealVar *p = new RooRealVar(pName.c_str(), pName.c_str(), -0.05 * i, -0.5, 0.0);
    RooRealVar *cp = new RooRealVar(cpName.c_str(), cpName.c_str(), 0.9, 0., 1.);
    
    // 设置第一个cp参数的初始值
    if (i == 0) {
      if (order == 1) {
        cp->setVal(1.0);
      } else {
        cp->setVal(0.1);
      }
    }
    
    formulaArgs.add(*p);
    formulaArgs.add(*cp);
    
    // 构建归一化的指数函数公式
    // 格式: cp * exp(p * x) / ∫ exp(p * x) dx (从xmin到xmax)
    if (i > 0) {
      formula += "+";
    }
    
    // 归一化的指数函数公式
    formula += Form("@%d*TMath::Exp(@%d*@0)/((TMath::Exp(@%d*%d)-TMath::Exp(@%d*%d))/@%d)",
      formulaArgs.getSize() - 1,  // cp
      formulaArgs.getSize() - 2,  // p
      formulaArgs.getSize() - 2,  // p (用于积分上限)
      xmax,
      formulaArgs.getSize() - 2,  // p (用于积分下限)
      xmin,
      formulaArgs.getSize() - 2   // p (分母)
    );
  }
  
  formula += ")";
  cout << "Exponential PDF formula: " << formula << endl;
  
  // 创建指数PDF
  RooGenericPdf *expPdf = new RooGenericPdf(
    Form("%s_exp%d", prefix.c_str(), order),
    Form("%s_exp%d", prefix.c_str(), order),
    formula.c_str(),
    formulaArgs
  );
  
  // 如果需要盲化处理
  if (doBlind) {
    string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
    RooGenericPdf *blindPdf = new RooGenericPdf(
      Form("%s_blind", prefix.c_str()),
      Form("%s_blind", prefix.c_str()),
      formula_SB.c_str(),
      RooArgList(*obs_var)
    );
    
    RooProdPdf *expPdf_SB = new RooProdPdf(
      Form("%s_SB", prefix.c_str()),
      "Blinded PDF (exclude 120-130)",
      RooArgList(*expPdf, *blindPdf)
    );
    
    return expPdf_SB;
  }
  
  return expPdf;
}


RooAbsPdf* PdfModelBuilder::getLaurentSeries(string prefix, int order){
 
  int nlower=int(ceil(order/2.));
  int nhigher=order-nlower;
  // first do 0th order
  RooArgList *pows = new RooArgList();
  RooArgList *plist = new RooArgList();
  string pname =  Form("%s_pow0",prefix.c_str());
  utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.))));
  pows->add(*utilities[pname]);

  // even terms
  for (int i=1; i<=nlower; i++){
    string name = Form("%s_l%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powl%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.-i))));
    pows->add(*utilities[pname]);
  }
  // odd terms
  for (int i=1; i<=nhigher; i++){
    string name = Form("%s_h%d",prefix.c_str(),i);
    params.insert(pair<string,RooRealVar*>(name, new RooRealVar(name.c_str(),name.c_str(),0.25/order,0.000001,0.999999)));
    plist->add(*params[name]);
    string pname =  Form("%s_powh%d",prefix.c_str(),i);
    utilities.insert(pair<string,RooAbsPdf*>(pname, new RooPower(pname.c_str(),pname.c_str(),*obs_var,RooConst(-4.+i))));
    pows->add(*utilities[pname]);
  }
  RooAddPdf *pdf = new RooAddPdf(prefix.c_str(),prefix.c_str(),*pows,*plist,true);
  return pdf;
  //bkgPdfs.insert(pair<string,RooAbsPdf*>(pdf->GetName(),pdf));
}


RooAbsPdf* PdfModelBuilder::getExpModGaussian(string prefix, bool doBlind){
  // Define parameters: mu, sigma and lambda
  RooRealVar *mu = new RooRealVar(Form("%s_mu",prefix.c_str()), Form("%s_mu",prefix.c_str()), 105., 95., 120.);
  RooRealVar *sigma = new RooRealVar(Form("%s_sigma",prefix.c_str()), Form("%s_sigma",prefix.c_str()), 3., 1., 10.); // Adjusted lower limit
  RooRealVar *lambda = new RooRealVar(Form("%s_lambda",prefix.c_str()), Form("%s_lambda",prefix.c_str()), 0.1, 0.001, 0.5); // Keep upper limit limited for stability
  
  // Add parameters to the params map
  params.insert(pair<string,RooRealVar*>(Form("%s_mu",prefix.c_str()), mu));
  params.insert(pair<string,RooRealVar*>(Form("%s_sigma",prefix.c_str()), sigma));
  params.insert(pair<string,RooRealVar*>(Form("%s_lambda",prefix.c_str()), lambda));
  
  RooExpModGauss *modgau = new RooExpModGauss(
      prefix.c_str(), "Exponentially Modified Gaussian", *obs_var, *mu,
      *sigma, *lambda);
    
  if (doBlind) {
    string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
    RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), formula_SB.c_str(), RooArgList(*obs_var));
    
    RooProdPdf* modgau_SB = new RooProdPdf(
        Form("%s_SB", prefix.c_str()),
        "Blinded PDF (exclude 120-130)",
        RooArgList(*modgau, *blindPdf)
    );
    return modgau_SB;
  }
 
  return modgau;
}


RooAbsPdf* PdfModelBuilder::getAsymGenGaussian(string prefix, bool doBlind){
  // Define parameters: beta, alpha and mu
  // This function implements the formula: beta/(2*alpha*Gamma(1/beta)) * exp(-((|x-mu|/alpha)^beta))
  // where:
  // - beta is the shape parameter (controls the peakedness)
  // - alpha is the scale parameter (controls the width)
  // - mu is the location parameter (controls the center position)
  
  RooRealVar *beta = new RooRealVar(Form("%s_beta", prefix.c_str()), Form("%s_beta", prefix.c_str()), 2.0, 0.5, 5.0);
  RooRealVar *alpha = new RooRealVar(Form("%s_alpha", prefix.c_str()), Form("%s_alpha", prefix.c_str()), 1.0, 0.1, 10.0);
  RooRealVar *mu = new RooRealVar(Form("%s_mu", prefix.c_str()), Form("%s_mu", prefix.c_str()), 105.0, 95.0, 125.0);
  
  // Add parameters to the params map
  params.insert(pair<string,RooRealVar*>(Form("%s_beta", prefix.c_str()), beta));
  params.insert(pair<string,RooRealVar*>(Form("%s_alpha", prefix.c_str()), alpha));
  params.insert(pair<string,RooRealVar*>(Form("%s_mu", prefix.c_str()), mu));
  
  // Create a list of formula variables
  RooArgList formulaVars;
  formulaVars.add(*obs_var);    // x
  formulaVars.add(*mu);         // mu
  formulaVars.add(*alpha);      // alpha
  formulaVars.add(*beta);       // beta
  
  // Formula implementation of: beta/(2*alpha*Gamma(1/beta)) * exp(-((|x-mu|/alpha)^beta))
  // We use TMath::Gamma(1/beta) for the gamma function
  // We need to ensure numerical stability:
  // 1. Prevent division by zero in gamma function
  // 2. Prevent overflow in the exponent
  // 3. Ensure the absolute value is properly handled
  string formula = "@3/(2*@2*TMath::Gamma(1.0/@3)) * exp(-pow(abs(@0-@1)/@2, @3))";
  
  // Create the generalized Gaussian PDF with the specified formula
  RooGenericPdf *agg = new RooGenericPdf(prefix.c_str(), "Generalized Gaussian Distribution", formula.c_str(), formulaVars);
      
  if (doBlind) {
    string formula_SB = "((@0 < 120)? 1:((@0 > 130)? 1:0))";
    RooGenericPdf *blindPdf = new RooGenericPdf(Form("%s_blind", prefix.c_str()), Form("%s_blind", prefix.c_str()), formula_SB.c_str(), RooArgList(*obs_var));
    
    RooProdPdf* agg_SB = new RooProdPdf(
        Form("%s_SB", prefix.c_str()),
        "Blinded PDF (exclude 120-130)",
        RooArgList(*agg, *blindPdf)
    );
    return agg_SB;
  }

  return agg;
}
