#include "../interface/GofBinIntegrals.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooNumIntConfig.h"
#include <cassert>
#include <iostream>

int main() {
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-12);
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-12);
  RooRealVar x("test_x", "x", .37, 0., 1.);
  RooRealVar a("test_a", "a", -.5, -5., 5.);
  RooExponential exponential("test_exponential", "exponential", x, a);
  int checks = 0;
  for (const auto& domain : {std::pair<double,double>(0.,1.), {-2.,3.}}) {
    x.setRange(domain.first, domain.second);
    for (double slope : {-2., -.01, 1.}) {
      a.setVal(slope);
      const double lo = domain.first + .2*(domain.second-domain.first);
      const double hi = domain.first + .45*(domain.second-domain.first);
      const double expected = (std::exp(slope*hi)-std::exp(slope*lo)) /
          (std::exp(slope*domain.second)-std::exp(slope*domain.first));
      for (const char* mode : {"raw", "normalized", "auto"}) {
        x.setVal(.37);
        {
          GofBinIntegrals calc(exponential, x, mode);
          assert(calc.usesRaw() == (std::string(mode)=="raw"));
          assert(std::fabs(calc.integral("test_bin",lo,hi)-expected)<1e-10);
          ++checks;
        }
        assert(x.getVal()==.37);
      }
    }
  }
  RooFormulaVar turnon("test_turnon", "(1+erf(@0))/2", RooArgList(x));
  RooEffProd direct("env_pdf_cat0_expbernerf2", "direct", exponential, turnon);
  assert(GofBinIntegrals::supportsAutomaticRaw(direct));
  direct.SetName("env_pdf_cat0_expcheberf2");
  assert(GofBinIntegrals::supportsAutomaticRaw(direct));
  direct.SetName("env_pdf_cat0_expcheberf2_SB");
  assert(!GofBinIntegrals::supportsAutomaticRaw(direct));
  direct.SetName("env_pdf_cat0_powerf1");
  assert(!GofBinIntegrals::supportsAutomaticRaw(direct));
  RooGenericPdf zero("test_zero", "zero", "0*@0", RooArgList(x));
  GofBinIntegrals fallback(zero,x,"raw");
  assert(!fallback.usesRaw() && fallback.fallbacks()==1);
  std::cout << "GOF_INTEGRALS_CHECK_PASS analytic=" << checks
            << " dispatch=4 state_restore=18 fallback=1" << std::endl;
}
