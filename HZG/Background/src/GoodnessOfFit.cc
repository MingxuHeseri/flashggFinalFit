#include "../interface/GoodnessOfFit.h"
#include "../interface/GofBinIntegrals.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <limits>
#include <vector>

#include "../interface/FitUtils.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TArrow.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TIterator.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

using namespace std;
using namespace RooFit;

namespace {
int countFloatingShapeParameters(RooAbsPdf *pdf, RooDataSet *data) {
  if (!pdf || !data) return 0;
  RooArgSet *params = pdf->getParameters(*data);
  if (!params) return 0;

  int nFloat = 0;
  for (RooAbsArg *obj : *params) {
    RooRealVar *var = dynamic_cast<RooRealVar*>(obj);
    if (!var) continue;
    if (var->isConstant()) continue;
    ++nFloat;
  }
  delete params;
  return nFloat;
}

double getBinnedGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data,
                              int nBinsForMass, const BackgroundFitOptions& options,
                              bool sidebandOnly, double *chi2Value, int *ndofValue) {
  static int histCounter = 0;
  const int histId = histCounter++;
  const double mggLow = mass->getMin();
  const double mggHigh = mass->getMax();
  TH1 *dataHist = data->createHistogram(Form("gof_data_%d", histId), *mass,
                                        RooFit::Binning(nBinsForMass, mggLow, mggHigh));
  if (!dataHist) {
    std::cerr << "[ERROR] Failed to create GOF data histogram" << std::endl;
    if (chi2Value) *chi2Value = 0.;
    if (ndofValue) *ndofValue = 0;
    return 0.;
  }

  const auto integrationStart = std::chrono::steady_clock::now();
  GofBinIntegrals binIntegrator(*mpdf, *mass, options.gofIntegralMode);
  std::vector<double> pdfBinIntegrals(dataHist->GetNbinsX() + 1, 0.);
  double dataYield = 0.;
  double pdfIntegral = 0.;
  int nUsedBins = 0;
  bool weightedData = false;

  for (int i = 1; i <= dataHist->GetNbinsX(); ++i) {
    const double binLow = dataHist->GetXaxis()->GetBinLowEdge(i);
    const double binHigh = dataHist->GetXaxis()->GetBinUpEdge(i);
    if (sidebandOnly && !((binHigh <= options.blindLow) || (binLow >= options.blindHigh))) continue;

    const TString rangeName = Form("gof_bin_%d_%d", histId, i);
    const double integral = binIntegrator.integral(rangeName.Data(), binLow, binHigh);

    if (integral <= 0.) continue;
    pdfBinIntegrals[i] = integral;
    pdfIntegral += integral;
    dataYield += dataHist->GetBinContent(i);
    const double observed = dataHist->GetBinContent(i);
    const double error2 = dataHist->GetBinError(i) * dataHist->GetBinError(i);
    if (observed > 0. && std::fabs(error2 - observed) > std::max(1.e-6, 1.e-6 * observed)) weightedData = true;
    ++nUsedBins;
  }

  std::cout << "[GOF_INTEGRATION] pdf=" << mpdf->GetName()
            << " requested=" << options.gofIntegralMode
            << " backend=" << (binIntegrator.usesRaw() ? "raw_shared_norm" : "normalized")
            << " fallback_bins=" << binIntegrator.fallbacks()
            << " seconds=" << std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - integrationStart).count() << std::endl;
  const int np = countFloatingShapeParameters(mpdf, data);
  int ndof = nUsedBins - np - 1;
  double pearsonChi2 = 0.;
  double poissonDeviance = 0.;
  double minimumExpected = std::numeric_limits<double>::infinity();
  double prob = 0.;
  double pearsonProb = 0.;
  double poissonProb = std::numeric_limits<double>::quiet_NaN();

  if (dataYield > 0. && pdfIntegral > 0. && ndof > 0) {
    for (int i = 1; i <= dataHist->GetNbinsX(); ++i) {
      if (pdfBinIntegrals[i] <= 0.) continue;
      const double observed = dataHist->GetBinContent(i);
      const double expected = dataYield * pdfBinIntegrals[i] / pdfIntegral;
      if (expected <= 0.) continue;
      minimumExpected = std::min(minimumExpected, expected);
      const double dataError = dataHist->GetBinError(i);
      const double variance = weightedData && dataError > 0.
                                  ? dataError * dataError : expected;
      if (variance > 0.) {
        pearsonChi2 += (observed - expected) * (observed - expected) / variance;
      }
      if (!weightedData) {
        double contribution = 2. * expected;
        if (observed > 0.) {
          contribution = 2. * (expected - observed + observed * std::log(observed / expected));
        }
        poissonDeviance += std::max(0., contribution);
      }
    }
    pearsonProb = TMath::Prob(pearsonChi2, ndof);
    if (weightedData) {
      prob = pearsonProb;
    } else {
      poissonProb = TMath::Prob(poissonDeviance, ndof);
      prob = std::min(pearsonProb, poissonProb);
    }
  } else {
    std::cerr << "[WARNING] GOF has invalid inputs: data yield = " << dataYield
              << ", pdf integral = " << pdfIntegral
              << ", used bins = " << nUsedBins
              << ", parameters = " << np << std::endl;
  }

  std::cout << "[INFO] " << (sidebandOnly ? "Sideband " : "") << "GOF used bins = " << nUsedBins
            << ", data yield = " << dataYield
            << ", pdf integral = " << pdfIntegral
            << ", floating shape parameters = " << np
            << ", fitted normalization constraint = 1"
            << ", ndof = " << ndof
            << ", statistic = " << (weightedData ? "Pearson chi2" : "conservative min(Pearson,Poisson)") << std::endl;
  std::cout << "[GOF_METRICS] pdf=" << mpdf->GetName()
            << " bins=" << nUsedBins << " ndof=" << ndof
            << " pearson_chi2=" << pearsonChi2
            << " pearson_prob=" << pearsonProb
            << " poisson_deviance=" << (weightedData ? -1. : poissonDeviance)
            << " poisson_prob=" << (weightedData ? -1. : poissonProb)
            << " conservative_prob=" << prob
            << " minimum_expected=" << minimumExpected << std::endl;
  if (chi2Value) *chi2Value = weightedData ? pearsonChi2 : std::max(pearsonChi2, poissonDeviance);
  if (ndofValue) *ndofValue = ndof;
  delete dataHist;
  return prob;
}
}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, std::string name, BackgroundFitOptions& options, bool sidebandOnly, double *chi2Value, int *ndofValue){
  const int nBinsForMass = std::max(
      1, static_cast<int>(std::round(
             (mass->getMax() - mass->getMin()) * options.gofBinsPerGeV)));

  // Routine to calculate the goodness of fit. 
  name+="_gofTest.pdf";
  const int np = countFloatingShapeParameters(mpdf, data);
  std::cout << "[INFO] Calculating " << (sidebandOnly ? "sideband " : "")
            << "GOF for pdf " << mpdf->GetName()
            << ", using " << np << " fitted parameters" <<std::endl;

  print_params(mpdf);

  double prob = getBinnedGoodnessOfFit(mass, mpdf, data, nBinsForMass, options, sidebandOnly, chi2Value, ndofValue);
  std::cout << "[INFO] Chi2 in " << (sidebandOnly ? "Sideband" : "Observed")
            << " =  " << (chi2Value ? *chi2Value : 0.) << std::endl;
  std::cout << "[INFO] p-value  =  " << prob << std::endl;
  return prob;

}

/**
 * Calculates the Kolmogorov-Smirnov test probability for a PDF fitted to data
 * @param mass - The mass variable
 * @param pdf - The PDF to test
 * @param data - The dataset to test against
 * @param name - Base name for output files
 * @return KS test probability
 */
double getKSProb(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, const BackgroundFitOptions& options) {
  // Calculate bin density based on data entries per GeV
  int mgg_low = mass->getMin(), mgg_high = mass->getMax();
  std::cout << "[INFO] KS test mass range: " << mgg_low << " - " << mgg_high << std::endl;

  // Determine bin density based on event statistics
  const double binForGeV = options.gofBinsPerGeV;
  const int nBin = std::max(1, static_cast<int>(std::round(
      (mgg_high - mgg_low) * binForGeV)));

  // Create histograms for data and PDF
  TH1F *dataHist = NULL;
  
  // Handle blinded or unblinded data
  if (options.blind) {
    // Create separate datasets for regions outside blind region
    RooDataSet *lowSideband = (RooDataSet*)data->reduce(Form("CMS_hzg_mass < %f", options.blindLow));
    RooDataSet *highSideband = (RooDataSet*)data->reduce(Form("CMS_hzg_mass > %f", options.blindHigh));
    
    // Convert to histograms
    TH1F *dataHistLow = (TH1F*)lowSideband->createHistogram("dataHistLow", *mass, 
                         RooFit::Binning(std::max(1, static_cast<int>(std::round((options.blindLow-mgg_low) * binForGeV))), mgg_low, options.blindLow));
    TH1F *dataHistHigh = (TH1F*)highSideband->createHistogram("dataHistHigh", *mass, 
                          RooFit::Binning(std::max(1, static_cast<int>(std::round((mgg_high-options.blindHigh) * binForGeV))), options.blindHigh, mgg_high));
    
    // Combine into one histogram
    dataHist = new TH1F("dataHist", "Data Histogram", nBin, mgg_low, mgg_high);
    
    // Fill with data from sidebands
    for (int i = 1; i <= dataHistLow->GetNbinsX(); i++) {
      dataHist->SetBinContent(i, dataHistLow->GetBinContent(i));
      dataHist->SetBinError(i, dataHistLow->GetBinError(i));
    }
    
    int offsetBin = static_cast<int>(std::round((options.blindHigh - mgg_low) * binForGeV));
    for (int i = 1; i <= dataHistHigh->GetNbinsX(); i++) {
      dataHist->SetBinContent(i + offsetBin, dataHistHigh->GetBinContent(i));
      dataHist->SetBinError(i + offsetBin, dataHistHigh->GetBinError(i));
    }
    
    delete dataHistLow;
    delete dataHistHigh;
    delete lowSideband;
    delete highSideband;
  } else {
    // Use all data if not blinded
    dataHist = (TH1F*)data->createHistogram("dataHist", *mass, RooFit::Binning(nBin, mgg_low, mgg_high));
  }
  
  // Generate histogram from PDF
  TH1F *pdfHist = (TH1F*)pdf->createHistogram("pdfHist", *mass, RooFit::Binning(nBin, mgg_low, mgg_high));
  // Set bin errors to zero
  if (pdfHist) {
    for (int i = 0; i <= pdfHist->GetNbinsX() + 1; ++i) { // Include underflow/overflow bins just in case
      pdfHist->SetBinError(i, 0.0);
    }
  }
  
  // Set signal region to zero if blinded
  if (options.blind) {
    int lowBin = pdfHist->FindBin(options.blindLow);
    int highBin = pdfHist->FindBin(options.blindHigh);
    for (int i = lowBin; i < highBin; i++) {
      pdfHist->SetBinContent(i, 0);
      dataHist->SetBinContent(i, 0);
    }
  }

  // Normalize histograms
  if (pdfHist->Integral() > 0) pdfHist->Scale(1.0/pdfHist->Integral());
  if (dataHist->Integral() > 0) dataHist->Scale(1.0/dataHist->Integral());
  
  // Get cumulative distributions
  TH1F *dataCum = (TH1F*)dataHist->GetCumulative();
  TH1F *pdfCum = (TH1F*)pdfHist->GetCumulative();

  // Normalize cumulative distributions
  if (dataCum->Integral() > 0) dataCum->Scale(1.0/dataCum->GetBinContent(dataCum->GetNbinsX()));
  if (pdfCum->Integral() > 0) pdfCum->Scale(1.0/pdfCum->GetBinContent(pdfCum->GetNbinsX()));
  
  // Calculate KS test probability using the histograms directly
  double ksProb = dataHist->KolmogorovTest(pdfHist); // Use TH1::KolmogorovTest

  // Create plot
  TCanvas *canKS = new TCanvas("canKS", "canKS", 800, 800);
  
  // Define uniform text size
  Float_t axisLabelSize = 0.06;   // Axis label size
  Float_t axisTitleSize = 0.07;   // Axis title size
  Float_t axisTitleOffset = 0.9;  // Axis title offset
  Float_t legendTextSize = 0.06;  // Legend text size
  Float_t latexTextSize = 0.06;   // Other text annotation size
  
  // Create pad layout
  double padHeight = 0.4;
  TPad *pad1 = new TPad("pad1","pad1",0,padHeight,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,padHeight);
  double ratio = (1-padHeight)/padHeight;
  pad1->SetBottomMargin(0.001);
  pad1->SetLeftMargin(0.12);
  pad2->SetTopMargin(0.001);
  pad2->SetLeftMargin(0.12);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  
  // Top pad: draw distributions
  pad1->cd();
  dataHist->SetLineColor(kBlack);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerColor(kBlack);
  dataHist->SetTitle("");
  dataHist->GetXaxis()->SetTitle(""); // Remove X axis title from top plot
  dataHist->GetYaxis()->SetTitle("Normalized Events");
  
  // Set axis text size for top pad
  dataHist->GetYaxis()->SetTitleSize(axisTitleSize);
  dataHist->GetYaxis()->SetLabelSize(axisLabelSize);
  dataHist->GetYaxis()->SetTitleOffset(axisTitleOffset);
  dataHist->GetXaxis()->SetLabelSize(0); // Hide X axis labels on top pad
  
  gPad->SetLeftMargin(0.12);
  gStyle->SetOptTitle(0);

  // Ensure Y axis uses appropriate format
  double maxY = dataHist->GetMaximum();
  if (maxY > 0 && maxY < 0.1) {
    dataHist->GetYaxis()->SetMaxDigits(2);
    dataHist->GetYaxis()->CenterTitle(true);
  }
  
  dataHist->Draw("E");
  
  pdfHist->SetLineColor(kBlue);
  pdfHist->SetLineWidth(2);
  pdfHist->Draw("HIST SAME");
  
  // Add legend with uniform text size
  TLegend *leg = new TLegend(0.6, 0.7, 0.89, 0.89);
  leg->SetTextSize(legendTextSize);
  leg->AddEntry(dataHist, "Data", "EP");
  leg->AddEntry(pdfHist, "PDF", "L");
  leg->AddEntry((TObject*)0, Form("KS Prob = %.4f", ksProb), "");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  
  // Add blind region (if needed)
  if (options.blind) {
    TBox *blindBox = new TBox(options.blindLow, 0, options.blindHigh, dataHist->GetMaximum()*1.1);
    blindBox->SetFillColorAlpha(kGray, 1);
    blindBox->SetFillStyle(3004);
    blindBox->Draw("same");
    
    TLegend *blindLeg = new TLegend(0.6, 0.6, 0.89, 0.7);
    blindLeg->SetTextSize(legendTextSize);
    blindLeg->SetBorderSize(0);
    blindLeg->SetFillStyle(0);
    blindLeg->AddEntry(blindBox, "Blinded Region", "f");
    blindLeg->Draw();
  }
  
  // Bottom pad: draw cumulative distributions
  pad2->cd();
  dataCum->SetLineColor(kBlack);
  dataCum->SetMarkerStyle(20);
  dataCum->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  dataCum->GetYaxis()->SetTitle("Cum. Prob.");
  dataCum->GetYaxis()->SetRangeUser(0, 1.1);
  
  // Set axis text size for bottom pad, ensure consistency with top pad
  dataCum->GetXaxis()->SetTitleSize(axisTitleSize*ratio);
  dataCum->GetYaxis()->SetTitleSize(axisTitleSize*ratio);
  dataCum->GetXaxis()->SetLabelSize(axisLabelSize*ratio);
  dataCum->GetYaxis()->SetLabelSize(axisLabelSize*ratio);
  dataCum->GetYaxis()->SetTitleOffset(axisTitleOffset*ratio);
  dataCum->GetXaxis()->SetTitleOffset(1.0);
  
  dataCum->Draw("E");
  
  pdfCum->SetLineColor(kBlue);
  pdfCum->SetLineWidth(2);
  pdfCum->Draw("HIST SAME");
  
  // Draw maximum KS distance
  double maxDiff = 0;
  int maxBin = 0;
  for (int i = 1; i <= dataCum->GetNbinsX(); i++) {
    double diff = fabs(dataCum->GetBinContent(i) - pdfCum->GetBinContent(i));
    if (diff > maxDiff) {
      maxDiff = diff;
      maxBin = i;
    }
  }
  
  // Add arrow to show KS distance
  double xAtMax = dataCum->GetBinCenter(maxBin);
  double y1 = dataCum->GetBinContent(maxBin);
  double y2 = pdfCum->GetBinContent(maxBin);
  
  TArrow *arrow = new TArrow(xAtMax, y1, xAtMax, y2, 0.02, "|>");
  arrow->SetLineColor(kRed);
  arrow->SetLineWidth(2);
  arrow->Draw();
  
  // Add label with uniform text size
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(latexTextSize*ratio);
  lat->DrawLatex(0.15, 0.85, Form("KS statistic = %.4f", maxDiff));
  
  // Add blind region to cumulative plot (if needed)
  if (options.blind) {
    TBox *blindBoxCum = new TBox(options.blindLow, 0, options.blindHigh, 1.1);
    blindBoxCum->SetFillColorAlpha(kGray, 1);
    blindBoxCum->SetFillStyle(3004);
    blindBoxCum->Draw("same");
  }
  
  // Save image
  canKS->SaveAs(Form("%s_KSTest.png", name.c_str()));
  
  // Clean up
  delete canKS;
  delete dataHist;
  delete pdfHist;
  delete dataCum;
  delete pdfCum;
  delete arrow;
  delete lat;
  delete leg;
  
  return ksProb;
}
