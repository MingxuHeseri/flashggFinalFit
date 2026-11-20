#include "../interface/FitUtils.h"
#include "../interface/SidebandChi2.h"
#include "RooBinning.h"
#include <algorithm>
#include <stdexcept>
#include <vector>

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "../interface/PdfFactory.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooMinimizer.h"
#include "TIterator.h"

using namespace std;
using namespace RooFit;

RooFitResult* Chi2Fit(RooAbsReal& chi2_, int printLevel_, double eps_, bool offSet_, int strategy_){
  // Your minimization settings (assuming these variables are defined)
  // int printLevel_ = -1; // Example value
  // double eps_ = 100;    // Example value
  // bool offSet_ = false;  // Example value
  // int strategy_ = 0;   // Example value
  // Create the minimizer instance using the chi2 variable calculated in sidebands
  RooMinimizer mini(chi2_);
  mini.setPrintLevel(printLevel_);
  mini.setEps(eps_);
  mini.setOffsetting(offSet_);
  mini.setStrategy(strategy_);
  // Perform the minimization
  mini.minimize("Minuit2", "migrad"); // Or "Minuit" if preferred/needed
  mini.hesse(); // Calculate uncertainties using HESSE
  // Optional: Get the fit result object
  RooFitResult *fitTest = mini.save("fitResult_Sidebands", "Fit Result from Sideband Chi2");
  // if (fitTest) {
  //     fitTest->Print("v"); // Print detailed fit result
  // }
  return fitTest;
}

void print_params(RooAbsPdf *pdf){
  RooArgSet* params = pdf->getParameters(RooArgSet());
  std::cout << "\n" << pdf->GetName() << std::endl;
  for (auto arg : *params) {
    RooRealVar* param = dynamic_cast<RooRealVar*>(arg);
    if (!param) continue;
    if (!param->isConstant()) {
      std::cout << param->GetName() << " = " << param->getVal() 
                << " (Range: [" << param->getMin() << ", " << param->getMax() << "])" << std::endl;
    }
    else {
      std::cout << param->GetName() << " = " << param->getVal() << std::endl;
    }
  }
  delete params;
}

namespace {
Bool_t useSumW2ErrorForFit(const BackgroundFitOptions& options) {
  return options.isData ? kFALSE : kTRUE;
}

std::unique_ptr<RooDataHist> binSidebandEvents(RooDataSet& data, const BackgroundFitOptions& options) {
  const auto& axis = static_cast<const RooRealVar&>(*data.get()->find("CMS_hzg_mass"));
  const auto& binning = axis.getBinning();
  std::vector<double> edges;
  for (int i=0; i<binning.numBins(); ++i) edges.push_back(binning.binLow(i));
  edges.push_back(binning.highBound());
  bool split = false;
  for (double edge : {options.mggLow, options.blindLow, options.blindHigh, options.mggHigh}) {
    if (edge < edges.front() || edge > edges.back())
      throw std::invalid_argument("Sidebands outside the data histogram domain");
    auto pos = std::lower_bound(edges.begin(), edges.end(), edge);
    if ((pos != edges.end() && std::abs(*pos-edge)<1e-9) ||
        (pos != edges.begin() && std::abs(*(pos-1)-edge)<1e-9)) continue;
    edges.insert(pos,edge);
    split = true;
  }
  if (!split && data.get()->size()==1) return std::unique_ptr<RooDataHist>(data.binnedClone());
  RooRealVar coordinate(axis.GetName(),axis.GetTitle(),edges.front(),edges.back());
  coordinate.setBinning(RooBinning(edges.size()-1,edges.data()));
  RooArgSet vars(coordinate);
  return std::make_unique<RooDataHist>("sidebandBins","Sideband events",vars,data);
}
}

/**
 * Perform fit and calculate chi2, save parameters if necessary
 * Supports both RooDataSet and RooDataHist types of data
 */
void runFit(RooAbsPdf *pdf, RooDataSet *data, RooDataSet *data_SB,
            double *NLL, int *stat_t, const BackgroundFitOptions& options,
            bool doBlind, const std::string& cat, const std::string& typePrefix,
            int order, int nBins, const RooAbsPdf* pdfToSave, bool chi2fit,
            int* covQualOut, double* edmOut, int* boundaryHitOut) {
  if (!doBlind) std::cout << "[INFO] Perform full range fit for RooDataSet" << std::endl;
  else std::cout << "[INFO] Perform sideband fit for RooDataSet" << std::endl;
  // Initialize variables
  int ntries = 0;
  int stat = 1;
  double minnll = 10e8;
  int lastCovQual = -999;
  double lastEdm = std::numeric_limits<double>::infinity();
  std::unique_ptr<RooArgSet> params_pdf(pdf->getParameters(*data));
  print_params(pdf);

  const char* sidebandRangeNames = "fitSidebandLow,fitSidebandHigh";
  auto* mass = dynamic_cast<RooRealVar*>(data->get()->find("CMS_hzg_mass"));
  if (doBlind) {
    std::unique_ptr<RooArgSet> obs(pdf->getObservables(*data_SB->get()));
    auto* pdfMass = dynamic_cast<RooRealVar*>(obs->find("CMS_hzg_mass"));
    auto* dataMass = dynamic_cast<RooRealVar*>(data_SB->get()->find("CMS_hzg_mass"));
    if (!mass || !pdfMass || !dataMass)
      throw std::invalid_argument("Sideband fit requires CMS_hzg_mass in PDF and data");
    for (auto* axis : {mass, pdfMass, dataMass}) {
      axis->setRange("fitSidebandLow", options.mggLow, options.blindLow);
      axis->setRange("fitSidebandHigh", options.blindHigh, options.mggHigh);
    }
  }
  std::unique_ptr<RooDataHist> sidebandDataHist;
  auto makeBinnedSidebandChi2 = [&]() {
    if (!sidebandDataHist) sidebandDataHist = binSidebandEvents(*data_SB, options);
    std::cout << "[INFO] Sideband chi2 mode=" << options.sidebandChi2Mode
              << " core=" << pdf->GetName() << std::endl;
    return SidebandChi2::make(*pdf,*sidebandDataHist,*mass,options.mggLow,options.blindLow,
                              options.blindHigh,options.mggHigh,options.sidebandChi2Mode=="integral");
  };

  // Perform fit
  while (stat != 0 && ntries < options.maxTries) {
    RooFitResult *fitTest = nullptr;
    
    if (doBlind && mass) {
      if (!chi2fit) {
        std::cout << "[INFO] Perform native two-range sideband fit for "
                  << pdf->GetName() << std::endl;
        fitTest = pdf->fitTo(*data_SB,
                            RooFit::Save(1),
                            RooFit::Range(sidebandRangeNames),
                            RooFit::Minimizer("Minuit2","migrad"),
                            RooFit::Strategy(options.minimizerStrategy),
                            RooFit::Offset(true),
                            RooFit::Hesse(true),
                            RooFit::MaxCalls(100000),
                            RooFit::SumW2Error(useSumW2ErrorForFit(options)),
                            RooFit::PrintLevel(-1));
      } else {
        std::cout << "[INFO] Using binned chi2 sideband fit for "
                  << pdf->GetName() << std::endl;
        auto chi2 = makeBinnedSidebandChi2();
        fitTest = Chi2Fit(*chi2, -1, options.minimizerEps, true,
                          options.minimizerStrategy);
      }
      
    } else {
      // Regular fit
      fitTest = pdf->fitTo(*data,
                          RooFit::Save(1),
                          RooFit::Minimizer("Minuit2","migrad"),
                          RooFit::Strategy(options.minimizerStrategy),
                          RooFit::Offset(true),
                          RooFit::Hesse(true),
                          RooFit::MaxCalls(100000),
                          RooFit::SumW2Error(useSumW2ErrorForFit(options)),
                          RooFit::PrintLevel(-1));
      // RooDataHist* data_hist = data->binnedClone();
      // RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data_hist, RooFit::Save(true), RooFit::SumW2Error(kTRUE), RooFit::Extended(true));
      // fitTest = Chi2Fit(chi2_);

    }
    
    if (fitTest) {
      stat = fitTest->status();
      minnll = fitTest->minNll();
      lastCovQual = fitTest->covQual();
      lastEdm = fitTest->edm();

      if (stat != 0 && ntries+1 < options.maxTries)
        params_pdf->assignValueOnly(fitTest->randomizePars());
      if (stat == 0 && nBins > 0 && !cat.empty() && mass) {
        print_params(pdf);
        const RooAbsPdf* targetPdf = pdfToSave ? pdfToSave : pdf;
        saveParamsToJSON(cat, typePrefix, order, targetPdf);
      }

      delete fitTest;
    } 
    else {
      stat = 1; // Mark as failed
    }
    
    ntries++;
  }

  // Report matched objectives once at the final parameter point.  The
  // likelihood is deliberately evaluated on the event-level sideband data,
  // while chi2 uses the cached explicitly binned representation.
  if (doBlind && mass && options.skipCrossObjectiveDiagnostics) {
    std::cout << "[CANDIDATE_OBJECTIVE_SKIPPED] category=" << cat
              << " family=" << typePrefix << " order=" << order
              << " reason=explicit_skip_cross_objective_diagnostics" << std::endl;
  }
  if (doBlind && mass && !options.skipCrossObjectiveDiagnostics) {
    std::unique_ptr<RooAbsReal> finalNll(pdf->createNLL(
        *data_SB, RooFit::Range(sidebandRangeNames), RooFit::Offset(false)));
    auto finalChi2 = makeBinnedSidebandChi2();
    // Evaluate before starting the marker line: RooFit may emit integration
    // diagnostics during getVal(), which must not corrupt machine parsing.
    const double finalNllValue = finalNll->getVal();
    const double finalChi2Value = finalChi2->getVal();
    const auto chi2Bins = finalChi2->bins().size();
    std::cout << "[CANDIDATE_OBJECTIVE] category=" << cat
              << " family=" << typePrefix << " order=" << order
              << " fit_objective=" << (chi2fit ? "binned_chi2" : "unbinned_nll")
              << " minimizer_value=" << minnll
              << " unbinned_nll=" << finalNllValue
              << " binned_chi2=" << finalChi2Value
              << " chi2_bins=" << chi2Bins
              << std::endl;
  }
  // Common unoffset NLL for tail-only order comparisons; legacy unchanged.
  if (doBlind && mass && !chi2fit &&
      (typePrefix == "exppolytail" || typePrefix == "powerlogtail" ||
       typePrefix == "laurenttail")) {
    std::unique_ptr<RooAbsReal> commonNll(pdf->createNLL(
        *data_SB, RooFit::Range(sidebandRangeNames), RooFit::Offset(false)));
    minnll = commonNll->getVal();
  }

  // Set output values
  *stat_t = stat;
  *NLL = minnll;

  int boundaryHit = 0;
  RooArgSet* fittedParams = pdf->getParameters(*data_SB);
  if (fittedParams) {
    for (auto argument : *fittedParams) {
      RooRealVar* parameter = dynamic_cast<RooRealVar*>(argument);
      if (!parameter) continue;
      if (parameter->isConstant()) continue;
      const double low = parameter->getMin();
      const double high = parameter->getMax();
      const double width = high - low;
      const double value = parameter->getVal();
      if (width > 0. &&
          (std::fabs(value - low) < 1.e-4 * width ||
           std::fabs(high - value) < 1.e-4 * width)) {
        boundaryHit = 1;
        std::cout << "[BOUNDARY_HIT] category=" << cat
                  << " parameter=" << parameter->GetName()
                  << " value=" << value
                  << " min=" << low << " max=" << high << std::endl;
      }
    }
    delete fittedParams;
  }
  if (covQualOut) *covQualOut = lastCovQual;
  if (edmOut) *edmOut = lastEdm;
  if (boundaryHitOut) *boundaryHitOut = boundaryHit;
}

// void runFit(RooAbsPdf *pdf, RooDataHist *data, double *NLL, int *stat_t, int MaxTries, bool doBlind=false, const std::string& cat = "", const std::string& typePrefix = "", int order = -1, int nBins = -1, const RooAbsPdf* pdfToSave = nullptr) {
//   std::cout << "[INFO] Perform sideband fit for RooDataHist" << std::endl;
//   // Initialize variables
//   int ntries = 0;
//   int stat = 1;
//   double minnll = 10e8;
//   RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
  
//   // chi2 calculation variables
//   double chi2_before = std::numeric_limits<double>::max();
//   RooRealVar* mass = nullptr;
//   std::string data_hist_name;
  
//   // Calculate chi2 before fit
//   if (nBins > 0 && !cat.empty()) {
//   mass = (RooRealVar*)data->get()->find("CMS_hgg_mass");
//   if (mass) {
//     RooPlot* frame_before = mass->frame();
//     data->plotOn(frame_before, Binning(nBins));
//     data_hist_name = frame_before->getHist()->GetName();
//     pdf->plotOn(frame_before, Name("pdf_before"));
    
//     int npars_before = pdf->getParameters(*data)->getSize();
//     chi2_before = frame_before->chiSquare("pdf_before", data_hist_name.c_str(), npars_before);
//     std::cout << "[INFO] Pre-fit chi2 = " << chi2_before << " for " << pdf->GetName() << std::endl;
    
//     delete frame_before;
//   } else {
//     nBins = -1; // Disable chi2 calculation if mass variable is not found
//   }
//   }

//   // Perform fit
//   while (stat != 0 && ntries < options.maxTries) {
//   RooFitResult *fitTest = nullptr;
  
//   if (doBlind && mass) {
//     // Save original range
//     double origMin = mass->getMin();
//     double origMax = mass->getMax();
    
//     // Set blind region range
//     mass->setRange("lowSideband", origMin, options.blindLow);
//     mass->setRange("highSideband", options.blindHigh, origMax);
    
//     // For RooDataHist, use Chi2 method to fit sideband regions
//     RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data, 
//             RooFit::Range("lowSideband,highSideband"), 
//             RooFit::Save(true), RooFit::SumW2Error(kTRUE));
//     fitTest = Chi2Fit(chi2_);
    
//     // Restore original range
//     mass->setRange(origMin, origMax);
//   } else {
//     // Regular fit - use Chi2 method for RooDataHist
//     RooChi2Var chi2_("chi2Var", "chi2Var", *pdf, *data, 
//             RooFit::Save(true), RooFit::SumW2Error(kTRUE));
//     fitTest = Chi2Fit(chi2_);
//   }
  
//   if (fitTest) {
//     stat = fitTest->status();
//     minnll = fitTest->minNll();

//     // Calculate chi2 after fit and save parameters if improved
//     if (stat == 0 && nBins > 0 && !cat.empty() && mass) {
//     RooPlot* frame_after = mass->frame();
//     data->plotOn(frame_after, Binning(nBins));
//     pdf->plotOn(frame_after, Name("pdf_after"));
    
//     int npars_after = pdf->getParameters(*data)->getSize();
//     double chi2_after = frame_after->chiSquare("pdf_after", data_hist_name.c_str(), npars_after);
    
//     // Save parameters if fit improved chi2
//     if (chi2_after < chi2_before) {
//       std::cout << "[INFO] Fit improved chi2: " << chi2_before << " -> " << chi2_after 
//           << " for " << typePrefix << order << std::endl;
//       saveParamsToJSON(cat, typePrefix, order, pdf);
//     }
    
//     delete frame_after;
//     }

//     // If fit failed, randomize parameters and try again
//     if (stat != 0) {
//     params_test->assignValueOnly(fitTest->randomizePars());
//     }
    
//     delete fitTest;
//   } else {
//     stat = 1; // Mark as failed
//   }
  
//   ntries++;
//   }
  
//   // Set output values
//   *stat_t = stat;
//   *NLL = minnll;
// }

// double getProbabilityFtest(double chi2, int ndof,RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, std::string name){
//   mgg_low = mass->getMin();
//   mgg_high = mass->getMax();
//   nBinsForMass = 4*(mgg_high-mgg_low);
 
//   double prob_asym = TMath::Prob(chi2,ndof);
//   if (!runFtestCheckWithToys) return prob_asym;

//   int ndata = data->sumEntries();
  
//   // fit the pdfs to the data and keep this fit Result (for randomizing)
//   RooFitResult *fitNullData = pdfNull->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
// 				,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME
//   RooFitResult *fitTestData = pdfTest->fitTo(*data,RooFit::Save(1),RooFit::Strategy(1)
// 				,RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE),RooFit::PrintLevel(-1)); //FIXME

//   // Ok we want to check the distribution in toys then 
//   // Step 1, cache the parameters of each pdf so as not to upset anything 
//   RooArgSet *params_null = pdfNull->getParameters((const RooArgSet*)(0));
//   RooArgSet preParams_null;
//   params_null->snapshot(preParams_null);
//   RooArgSet *params_test = pdfTest->getParameters((const RooArgSet*)(0));
//   RooArgSet preParams_test;
//   params_test->snapshot(preParams_test);
 
//   int ntoys =500;
//   TCanvas *can = new TCanvas();
//   can->SetLogy();
//   TH1F toyhist(Form("toys_fTest_%s.pdf",pdfNull->GetName()),";Chi2;",60,-2,10);
//   TH1I toyhistStatN(Form("Status_%s.pdf",pdfNull->GetName()),";FitStatus;",8,-4,4);
//   TH1I toyhistStatT(Form("Status_%s.pdf",pdfTest->GetName()),";FitStatus;",8,-4,4);

//   TGraph *gChi2 = new TGraph();
//   gChi2->SetLineColor(kGreen+2);
//   double w = toyhist.GetBinWidth(1);

//   int ipoint=0;

//   for (int b=0;b<toyhist.GetNbinsX();b++){
// 		double x = toyhist.GetBinCenter(b+1);
// 		if (x>0){
// 		  gChi2->SetPoint(ipoint,x,(ROOT::Math::chisquared_pdf(x,ndof)));
// 		  ipoint++;
// 		}
//   }
//   int npass =0; int nsuccesst =0;
//   mass->setBins(nBinsForMass);
//   for (int itoy = 0 ; itoy < ntoys ; itoy++){
//     params_null->assignValueOnly(preParams_null);
//     params_test->assignValueOnly(preParams_test);
//     RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass),ndata,0,1);

// 		int stat_n=1;
//         int stat_t=1;
// 		int ntries = 0;
// 		double nllNull,nllTest;
// 		// Iterate on the fit 
// 		int MaxTries = 3;
// 		while (stat_n!=0){
// 		  if (ntries>=MaxTries) break;
// 		  RooFitResult *fitNull = pdfNull->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
// 				,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
// 				//,RooFit::Optimize(0));

// 		  nllNull = fitNull->minNll();
//       stat_n = fitNull->status();
// 		  if (stat_n!=0) params_null->assignValueOnly(fitNullData->randomizePars());
// 		  ntries++; 
// 		}
		
// 		ntries = 0;
// 		while (stat_t!=0){
// 		  if (ntries>=MaxTries) break;
// 		  RooFitResult *fitTest = pdfTest->fitTo(*binnedtoy,RooFit::Save(1),RooFit::Strategy(1),RooFit::SumW2Error(kTRUE) //FIXME
// 				,RooFit::Minimizer("Minuit2","minimize"),RooFit::Minos(0),RooFit::Hesse(0),RooFit::PrintLevel(-1));
// 		  nllTest = fitTest->minNll();
//       stat_t = fitTest->status();
// 		  if (stat_t!=0) params_test->assignValueOnly(fitTestData->randomizePars()); 
// 		  ntries++; 
// 		}
       
// 		toyhistStatN.Fill(stat_n);
// 		toyhistStatT.Fill(stat_t);

//     if (stat_t !=0 || stat_n !=0) continue;
// 		nsuccesst++;
// 		double chi2_t = 2*(nllNull-nllTest);
// 		if (chi2_t >= chi2) npass++;
//     toyhist.Fill(chi2_t);
//   }

//   double prob=0;
//   if (nsuccesst!=0)  prob = (double)npass / nsuccesst;
//   toyhist.Scale(1./(w*toyhist.Integral()));
//   toyhist.Draw();
//   TArrow lData(chi2,toyhist.GetMaximum(),chi2,0);
//   lData.SetLineWidth(2);
//   lData.Draw();
//   gChi2->Draw("L");
//   TLatex *lat = new TLatex();
//   lat->SetNDC();
//   lat->SetTextFont(42);
//   lat->DrawLatex(0.1,0.91,Form("Prob (asymptotic) = %.4f (%.4f)",prob,prob_asym));
//   can->SaveAs(name.c_str());

//   TCanvas *stas =new TCanvas();
//   toyhistStatN.SetLineColor(2);
//   toyhistStatT.SetLineColor(1); 
//   TLegend *leg = new TLegend(0.2,0.6,0.4,0.87); leg->SetFillColor(0);
//   leg->SetTextFont(42);
//   leg->AddEntry(&toyhistStatN,"Null Hyp","L");
//   leg->AddEntry(&toyhistStatT,"Test Hyp","L");
//   toyhistStatN.Draw();
//   toyhistStatT.Draw("same");
//   leg->Draw();
//   stas->SaveAs(Form("%s_fitstatus.pdf",name.c_str()));
//   //reassign params
//   params_null->assignValueOnly(preParams_null);
//   params_test->assignValueOnly(preParams_test);

//   delete can; delete stas;
//   delete gChi2;
//   delete leg;
//   delete lat;

//   // Still return the asymptotic prob (usually its close to the toys one)
//   return prob_asym;

// }
