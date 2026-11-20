#ifndef FLASHGGFINALFIT_HZG_GOF_BIN_INTEGRALS_H
#define FLASHGGFINALFIT_HZG_GOF_BIN_INTEGRALS_H

#include <cmath>
#include <memory>
#include <string>
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooRealVar.h"

// One fixed-parameter GOF evaluation owns one normalization. No values survive
// across candidates, parameter changes, observable domains or separate calls.
// The RooFit integrator and its precision settings are deliberately unchanged.
class GofBinIntegrals {
 public:
  static bool supportsAutomaticRaw(const RooAbsPdf& pdf) {
    const std::string kind = pdf.ClassName();
    if (kind != "RooEffProd" && kind != "RooExpPolyErfPdf" && kind != "RooPowerErfPdf") return false;
    const std::string name = pdf.GetName();
    if (kind == "RooPowerErfPdf") {
      if (pdf.normRange()) return false;
      // Factory uses component count for PowerLaw and order for Laurent.
      for (const char* suffix : {"_powerf1", "_powerf2", "_powerf3", "_lauerf2",
                                 "_lauerf3", "_lauerf4", "_lauerf5", "_lauerf6"}) {
        if (name.ends_with(suffix)) return true;
      }
      return false;
    }
    // Only these two families have passed the production-snapshot precision
    // replay. Other direct-erf families retain their original integral path.
    for (const char* token : {"_expbernerf", "_expcheberf"}) {
      const auto pos = name.rfind(token);
      if (pos == std::string::npos) continue;
      const std::string order = name.substr(pos + std::string(token).size());
      if (!order.empty() && order.find_first_not_of("0123456789") == std::string::npos)
        return true;
    }
    return false;  // In particular, do not auto-select masked or FFT PDFs.
  }

  GofBinIntegrals(RooAbsPdf& pdf, RooRealVar& mass, const std::string& mode)
      : pdf_(pdf), mass_(mass), obs_(mass), savedMass_(mass.getVal()),
        raw_(mode == "raw" || (mode == "auto" && supportsAutomaticRaw(pdf))) {
    if (raw_) {
      std::unique_ptr<RooAbsReal> full(pdf_.createIntegral(obs_));
      norm_ = full ? full->getVal() : 0.;
      if (!(std::isfinite(norm_) && norm_ > 0.)) {
        raw_ = false;
        ++fallbacks_;
      }
    }
  }
  ~GofBinIntegrals() { mass_.setVal(savedMass_); }
  GofBinIntegrals(const GofBinIntegrals&) = delete;
  GofBinIntegrals& operator=(const GofBinIntegrals&) = delete;

  double integral(const char* rangeName, double low, double high) {
    mass_.setRange(rangeName, low, high);
    if (raw_) {
      std::unique_ptr<RooAbsReal> bin(pdf_.createIntegral(obs_, RooFit::Range(rangeName)));
      const double value = bin ? bin->getVal() / norm_ : -1.;
      if (std::isfinite(value) && value >= 0.) return value;
      ++fallbacks_;
    }
    std::unique_ptr<RooAbsReal> bin(pdf_.createIntegral(
        obs_, RooFit::NormSet(obs_), RooFit::Range(rangeName)));
    return bin ? bin->getVal() : 0.;
  }
  bool usesRaw() const { return raw_; }
  int fallbacks() const { return fallbacks_; }

 private:
  RooAbsPdf& pdf_;
  RooRealVar& mass_;
  RooArgSet obs_;
  double savedMass_;
  bool raw_;
  double norm_ = 1.;
  int fallbacks_ = 0;
};
#endif
