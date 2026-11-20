#ifndef FLASHGGFINALFIT_HZG_SIDEBAND_CHI2_H
#define FLASHGGFINALFIT_HZG_SIDEBAND_CHI2_H

#include "RooAbsReal.h"
#include "RooListProxy.h"
#include "RooRealProxy.h"
#include <memory>
#include <vector>

class RooAbsPdf;
class RooDataHist;
class RooRealVar;

// Shape-only chi-square. The two sidebands share a single normalization.
class SidebandChi2 final : public RooAbsReal {
public:
  struct Bin {
    double center, width, observed, errorLow, errorHigh;
  };
  static std::unique_ptr<SidebandChi2> make(const RooAbsPdf& core, RooDataHist& data,
      const RooRealVar& observable, double low, double blindLow, double blindHigh,
      double high, bool integrateBins);

  SidebandChi2(const SidebandChi2&, const char* name = nullptr);
  TObject* clone(const char* name) const override { return new SidebandChi2(*this, name); }
  const std::vector<Bin>& bins() const { return bins_; }
  double expected(std::size_t i) const;

private:
  SidebandChi2(RooAbsReal& core, RooRealVar& mass);
  double evaluate() const override;
  double numerator(std::size_t i) const;
  double scale() const;

  RooRealProxy core_, mass_;
  RooListProxy sides_, integrals_;
  std::vector<Bin> bins_;
  bool normalized_ = false;
  double total_ = 0.;
};
#endif
