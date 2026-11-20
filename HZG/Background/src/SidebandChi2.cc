#include "../interface/SidebandChi2.h"
#include "RooAbsBinning.h"
#include "RooAbsPdf.h"
#include "RooCustomizer.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooNumIntConfig.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

SidebandChi2::SidebandChi2(RooAbsReal& core, RooRealVar& mass)
    : RooAbsReal("sidebandChi2", "Sideband shape chi-square"),
      core_("core","shape",this,core), mass_("mass","private bin coordinate",this,mass),
      sides_("sides","sideband integrals",this), integrals_("integrals","bin integrals",this) {}

SidebandChi2::SidebandChi2(const SidebandChi2& other, const char* name)
    : RooAbsReal(other,name), core_("core",this,other.core_), mass_("mass",this,other.mass_),
      sides_("sides",this,other.sides_), integrals_("integrals",this,other.integrals_),
      bins_(other.bins_), normalized_(other.normalized_), total_(other.total_) {}

std::unique_ptr<SidebandChi2> SidebandChi2::make(const RooAbsPdf& core, RooDataHist& data,
    const RooRealVar& observable, double low, double blindLow, double blindHigh,
    double high, bool integrateBins) {
  if (!(low < blindLow && blindLow < blindHigh && blindHigh < high))
    throw std::invalid_argument("Invalid sideband intervals");
  const auto* coords = data.get();
  const auto* axis = dynamic_cast<const RooRealVar*>(coords->find(observable.GetName()));
  if (coords->size()!=1 || !axis) throw std::invalid_argument("Expected a one-dimensional mass histogram");
  std::unique_ptr<RooArgSet> deps(core.getObservables(*coords));
  auto* source = dynamic_cast<RooRealVar*>(deps->find(observable.GetName()));
  if (!source || source->getMin()>low || source->getMax()<high)
    throw std::invalid_argument("Sidebands outside the PDF mass domain");

  // RooCustomizer shares all shape parameters. Only the quadrature coordinate
  // is private, so named ranges and bin evaluations cannot mutate a caller's mass.
  // RooRealVar's copy constructor shares named binnings with its source.
  auto mass = std::make_unique<RooRealVar>("sidebandChi2_mass",source->GetTitle(),
      source->getVal(),source->getMin(),source->getMax(),source->getUnit());
  mass->setBinning(source->getBinning());
  if (source->hasBinning("cache")) mass->setBinning(source->getBinning("cache"), "cache");
  mass->setConstant(true);
  RooCustomizer custom(core,"sidebandChi2");
  custom.replaceArg(*source,*mass);
  std::unique_ptr<RooAbsPdf> shape(static_cast<RooAbsPdf*>(custom.build()));
  if (std::string(core.ClassName())=="RooExpPolyErfPdf") {
    // The same denominator scales every bin. Keep the accuracy of this
    // objective independent of a coarse global absolute integration tolerance.
    auto* config=shape->specialIntegratorConfig(true);
    config->setEpsAbs(std::min(config->epsAbs(),1e-50));
    config->setEpsRel(std::min(config->epsRel(),1e-10));
  }
  auto out = std::unique_ptr<SidebandChi2>(new SidebandChi2(*shape,*mass));
  const std::string kind = core.ClassName();
  out->normalized_ = kind!="RooErfProductPdf" && kind!="RooPowerErfPdf" && kind!="RooExpPolyErfPdf";
  RooArgSet obs(*mass);
  auto addIntegral = [&](RooListProxy& list, const std::string& name, double a, double b) {
    mass->setRange(name.c_str(),a,b);
    std::unique_ptr<RooAbsReal> term(out->normalized_
        ? shape->createIntegral(obs,RooFit::NormSet(obs),RooFit::Range(name.c_str()))
        : shape->createIntegral(obs,RooFit::Range(name.c_str())));
    list.add(*term);
    out->addOwnedComponents(std::move(term));
  };
  const auto& grid = axis->getBinning();
  for (double edge : {low,blindLow,blindHigh,high}) {
    // A reduced histogram cannot locate events within a bin cut by the mask.
    const int i = grid.binNumber(edge);
    if (std::abs(grid.binLow(i)-edge)>1e-9 && std::abs(grid.binHigh(i)-edge)>1e-9)
      throw std::invalid_argument("Sideband limits must coincide with histogram edges");
  }
  addIntegral(out->sides_,"chi2_low",low,blindLow);
  addIntegral(out->sides_,"chi2_high",blindHigh,high);
  for (int row=0; row<data.numEntries(); ++row) {
    data.get(row);
    const double center=axis->getVal();
    if (!(center>low && center<blindLow) && !(center>blindHigh && center<high)) {
      if (data.weight()!=0.) throw std::invalid_argument("Nonzero weight outside the declared sidebands");
      continue;
    }
    const int i=grid.binNumber(center);
    double errLow,errHigh;
    data.weightError(errLow,errHigh,RooAbsData::Poisson);
    out->bins_.push_back({center,grid.binWidth(i),data.weight(),errLow,errHigh});
    if (integrateBins)
      addIntegral(out->integrals_,"chi2_bin_"+std::to_string(i),grid.binLow(i),grid.binHigh(i));
  }
  out->total_=data.sumEntries();
  out->addOwnedComponents(std::move(shape),std::move(mass));
  return out;
}

double SidebandChi2::scale() const {
  const double den=static_cast<const RooAbsReal&>(sides_[0]).getVal()
                 + static_cast<const RooAbsReal&>(sides_[1]).getVal();
  return den>0. && std::isfinite(den) ? total_/den : std::numeric_limits<double>::quiet_NaN();
}

double SidebandChi2::numerator(std::size_t i) const {
  if (integrals_.size()) return static_cast<const RooAbsReal&>(integrals_[i]).getVal();
  static_cast<RooRealVar&>(*mass_.absArg()).setVal(bins_[i].center);
  if (normalized_) {
    RooArgSet obs(*mass_.absArg());
    return core_.arg().getVal(&obs)*bins_[i].width;
  }
  return core_.arg().getVal()*bins_[i].width;
}

double SidebandChi2::expected(std::size_t i) const {
  if (i>=bins_.size()) throw std::out_of_range("Sideband bin index");
  const double factor=scale();
  return factor*numerator(i);
}

double SidebandChi2::evaluate() const {
  const double factor=scale();
  double sum=0.,carry=0.;
  for (std::size_t i=0; i<bins_.size(); ++i) {
    const auto& b=bins_[i];
    const double mu=factor*numerator(i);
    const double delta=mu-b.observed;
    const double error=delta>0. ? b.errorHigh : b.errorLow;
    if (error==0. && mu==0. && b.observed==0.) continue;
    if (!(mu>=0. && std::isfinite(mu) && error>0. && std::isfinite(error))) {
      logEvalError("Undefined sideband expectation or Poisson error");
      return std::numeric_limits<double>::infinity();
    }
    const double corrected=delta*delta/(error*error)-carry;
    const double next=sum+corrected;
    carry=(next-sum)-corrected;
    sum=next;
  }
  return sum;
}
