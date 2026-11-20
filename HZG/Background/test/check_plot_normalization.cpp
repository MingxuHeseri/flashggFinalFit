#include "../interface/PlotNormalization.h"
#include "RooCurve.h"
#include "RooDataSet.h"
#include "RooEffProd.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooNumIntConfig.h"
#include "TROOT.h"
#include <cassert>
#include <iostream>
int main() {
  gROOT->SetBatch(true);
  RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-12);
  RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-12);
  RooRealVar x("x","x",.5,0.,1.), a("a","a",-.5,-5.,5.);
  RooExponential exp("exp","exp",x,a);
  RooFormulaVar eff("eff","(1+erf(@0))/2",RooArgList(x));
  RooEffProd pdf("test_expbernerf2","test",exp,eff);
  int cases=0, points=0;
  double maxDifference=0.;
  for (const auto& domain : {std::pair<double,double>(0.,1.), {-2.,3.}}) {
    x.setRange(domain.first,domain.second);
    const double width=domain.second-domain.first;
    x.setRange("low",domain.first,domain.first+.4*width);
    x.setRange("high",domain.first+.6*width,domain.second);
    pdf.setStringAttribute("fitrange","low,high");
    for (double slope : {-2.,-.01,1.}) {
      a.setVal(slope);
      for (bool withData : {false,true}) {
        RooArgSet obs(x);
        RooDataSet data("data","data",obs);
        for(int i=0;i<100;++i){
          if(i>=40 && i<60)continue;
          x.setVal(domain.first+(i+.5)*width/100);data.add(obs);
        }
        std::unique_ptr<RooPlot> old(x.frame()), fast(x.frame());
        if(withData){data.plotOn(old.get());data.plotOn(fast.get());}
        pdf.plotOn(old.get(),RooFit::LineColor(kRed),RooFit::LineStyle(2));
        x.setVal(.5);
        assert(plotWithSharedNormalization(pdf,x,*fast,kRed,2));
        assert(x.getVal()==.5);
        assert(a.getVal()==slope);
        assert(old->numItems()==fast->numItems());
        for(int i=withData?1:0;i<old->numItems();++i){
          auto* one=dynamic_cast<RooCurve*>(old->getObject(i));
          auto* two=dynamic_cast<RooCurve*>(fast->getObject(i));
          assert(one && two && one->GetN()==two->GetN());
          assert(one->GetLineColor()==two->GetLineColor());
          assert(one->GetLineStyle()==two->GetLineStyle());
          for(int j=0;j<one->GetN();++j){
            assert(one->GetPointX(j)==two->GetPointX(j));
            const double err=std::abs(one->GetPointY(j)-two->GetPointY(j));
            assert(err<1e-9+1e-7*std::abs(one->GetPointY(j)));
            maxDifference=std::max(maxDifference,err);++points;
          }
        }
        ++cases;
      }
    }
  }
  std::unique_ptr<RooPlot> frame(x.frame());
  pdf.setStringAttribute("fitrange",nullptr);
  assert(!plotWithSharedNormalization(pdf,x,*frame,kBlue,1));
  pdf.setStringAttribute("fitrange","missing");
  assert(!plotWithSharedNormalization(pdf,x,*frame,kBlue,1));
  pdf.setStringAttribute("fitrange","low,low");
  assert(!plotWithSharedNormalization(pdf,x,*frame,kBlue,1));
  pdf.setStringAttribute("fitrange","low,high");pdf.setNormRange("low");
  assert(!plotWithSharedNormalization(pdf,x,*frame,kBlue,1));pdf.setNormRange(nullptr);
  pdf.SetName("test_powerf1");
  assert(!plotWithSharedNormalization(pdf,x,*frame,kBlue,1));
  assert(frame->numItems()==0);
  std::cout<<"PLOT_NORMALIZATION_CHECK_PASS cases="<<cases<<" points="<<points
           <<" max_abs_difference="<<maxDifference<<" fallback=5"<<std::endl;
}
