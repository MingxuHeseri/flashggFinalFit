#ifndef HZG_PLOT_NORMALIZATION_H
#define HZG_PLOT_NORMALIZATION_H
#include <cmath>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooLinkedList.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "GofBinIntegrals.h"

// Plot-only optimization: preserve the ordinary normalized PDF projection and
// RooCurve sampler, but avoid integrating that projection a second time.
// Each call computes its own normalization; nothing is cached across fits.
inline bool plotWithSharedNormalization(RooAbsPdf& pdf, RooRealVar& mass,
                                        RooPlot& frame, int color, int style) {
  if (!GofBinIntegrals::supportsAutomaticRaw(pdf) || pdf.normRange()) return false;
  const char* fitRange = pdf.getStringAttribute("fitrange");
  if (!fitRange || !*fitRange) return false;
  if (frame.getPlotVar() == nullptr ||
      std::string(frame.getPlotVar()->GetName()) != mass.GetName()) return false;
  const double savedMass = mass.getVal();
  struct Restore { RooRealVar& x; double value; ~Restore(){x.setVal(value);} } restore{mass,savedMass};
  std::vector<std::pair<double,double>> ranges;
  std::stringstream names(fitRange);std::string name;
  double events = 0.;
  while (std::getline(names,name,',')) {
    if (name.empty() || !mass.hasRange(name.c_str())) return false;
    const double lo=mass.getMin(name.c_str()), hi=mass.getMax(name.c_str());
    if (!(std::isfinite(lo) && std::isfinite(hi) && lo < hi) ||
        lo < mass.getMin() || hi > mass.getMax()) return false;
    for (const auto& r:ranges) if (lo < r.second && hi > r.first) return false;
    ranges.emplace_back(lo,hi);
    events += frame.getFitRangeNEvt(lo,hi);
  }
  if (ranges.empty()) return false;
  RooArgSet obs(mass);
  std::unique_ptr<RooAbsReal> full(pdf.createIntegral(obs));
  std::unique_ptr<RooAbsReal> selected(pdf.createIntegral(obs,RooFit::Range(fitRange)));
  const double total=full->getVal(), sideband=selected->getVal();
  if (!(std::isfinite(total) && total>0 && std::isfinite(sideband) && sideband>0)) return false;
  // RooAbsPdf::plotOn uses event scaling only when the frame has data.
  const double scale=(frame.getFitRangeNEvt() ? events : 1.) *
                     frame.getFitRangeBinW() * total/sideband;
  if (!(std::isfinite(scale) && scale>0)) return false;
  frame.updateNormVars(mass);
  for (const auto& r:ranges) {
    auto range=RooFit::Range(r.first,r.second,false);
    auto norm=RooFit::Normalization(scale,RooAbsReal::Raw);
    auto lineColor=RooFit::LineColor(color);
    auto lineStyle=RooFit::LineStyle(style);
    RooLinkedList commands;
    commands.Add(&range);commands.Add(&norm);
    commands.Add(&lineColor);commands.Add(&lineStyle);
    // Call the RooAbsReal list overload explicitly: the variadic overload
    // redispatches to RooAbsPdf and would re-insert the automatic NormRange.
    pdf.RooAbsReal::plotOn(&frame,commands);
  }
  return true;
}
#endif
