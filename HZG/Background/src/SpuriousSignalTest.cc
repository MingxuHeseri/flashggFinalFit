#include "../interface/SpuriousSignalTest.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

#include "../interface/FitUtils.h"
#include "../interface/PdfFactory.h"
#include "../interface/PlotUtils.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooMinimizer.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"

using namespace std;
using namespace RooFit;

namespace {

constexpr const char* kSstWorkspaceName = "CMS_hzg_workspace";
constexpr const char* kMassVarName = "CMS_hzg_mass";
constexpr double kSignalPdfYieldScale = 1000.;

RooWorkspace* loadWorkspace(TFile* file, const TString& fileName) {
  if (!file || file->IsZombie()) return nullptr;
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(file->Get(kSstWorkspaceName));
  if (!ws) {
    std::cerr << "[ERROR] Workspace " << kSstWorkspaceName << " not found in " << fileName << std::endl;
  }
  return ws;
}

RooAbsData* loadDataset(RooWorkspace* ws, const std::vector<TString>& names, TString& selectedName) {
  if (!ws) return nullptr;
  for (const auto& name : names) {
    RooAbsData* data = ws->data(name.Data());
    if (data) {
      selectedName = name;
      return data;
    }
  }
  std::cerr << "[ERROR] None of the requested datasets were found:";
  for (const auto& name : names) std::cerr << " " << name;
  std::cerr << std::endl;
  return nullptr;
}

double datasetYield(RooAbsData* data, double fitLow, double fitHigh, double blindLow, double blindHigh) {
  if (!data) return 0.;

  double yield = 0.;
  for (int i = 0; i < data->numEntries(); ++i) {
    const RooArgSet* row = data->get(i);
    if (!row) continue;
    const double mass = row->getRealValue(kMassVarName);
    if (mass < fitLow || mass > fitHigh) continue;
    if (mass < blindLow || mass > blindHigh) yield += data->weight();
  }
  return yield;
}

TH1F* datasetToHist(RooAbsData* data, const TString& name, int nBins, double fitLow, double fitHigh) {
  auto* hist = new TH1F(name, name, nBins, fitLow, fitHigh);
  hist->Sumw2();

  if (!data) return hist;
  for (int i = 0; i < data->numEntries(); ++i) {
    const RooArgSet* row = data->get(i);
    if (!row) continue;
    const double mass = row->getRealValue(kMassVarName);
    if (mass < fitLow || mass > fitHigh) continue;
    hist->Fill(mass, data->weight());
  }
  return hist;
}

void smoothSstTemplate(TH1* hist) {
  if (!hist) return;
  const double rawIntegral = hist->Integral();
  if (!std::isfinite(rawIntegral) || rawIntegral <= 0.) {
    throw std::runtime_error("Cannot smooth an SST template with non-positive integral");
  }
  std::vector<double> rawErrors(hist->GetNbinsX() + 2, 0.);
  for (int bin = 0; bin <= hist->GetNbinsX() + 1; ++bin) {
    rawErrors[bin] = hist->GetBinError(bin);
  }
  hist->Smooth(1);
  const double smoothedIntegral = hist->Integral();
  if (!std::isfinite(smoothedIntegral) || smoothedIntegral <= 0.) {
    throw std::runtime_error("SST smoothing produced a non-positive integral");
  }
  const double correction = rawIntegral / smoothedIntegral;
  for (int bin = 0; bin <= hist->GetNbinsX() + 1; ++bin) {
    const double content = hist->GetBinContent(bin) * correction;
    if (!std::isfinite(content) || content < 0.) {
      throw std::runtime_error("SST smoothing produced a negative or non-finite bin");
    }
    hist->SetBinContent(bin, content);
    hist->SetBinError(bin, rawErrors[bin]);
  }
}

RooFitResult* fitSstNllTemplate(RooAbsPdf& pdf, RooAbsData& data,
                                bool extended, const TString& label) {
  std::unique_ptr<RooAbsReal> nll(
      pdf.createNLL(data, RooFit::Extended(extended)));
  if (!nll) return nullptr;
  RooMinimizer minimizer(*nll);
  minimizer.setPrintLevel(-1);
  minimizer.setStrategy(2);
  minimizer.setEps(1.e-4);
  minimizer.setOffsetting(true);
  minimizer.minimize("Minuit2", "migrad");
  minimizer.hesse();
  return minimizer.save(Form("fitResult_%s", label.Data()),
                        Form("SST Poisson fit %s", label.Data()));
}

RooFitResult* fitSstChi2Template(RooAbsPdf& pdf, RooDataHist& data,
                                 bool extended, const TString& label) {
  (void)label;
  std::unique_ptr<RooAbsReal> chi2(pdf.createChi2(
      data, RooFit::DataError(RooAbsData::SumW2),
      RooFit::Extended(extended)));
  if (!chi2) return nullptr;
  return Chi2Fit(*chi2, -1, 1.e-4, true, 2);
}

bool signalYieldFromPdfWorkspace(RooWorkspace* ws, const TString& cat,
                                 double& yield, TString& source) {
  if (!ws) return false;
  const TString eleName = Form("dcb_all_all_%s_ele_13TeV_norm", cat.Data());
  const TString muName = Form("dcb_all_all_%s_mu_13TeV_norm", cat.Data());
  RooAbsReal* eleNorm = ws->function(eleName.Data());
  RooAbsReal* muNorm = ws->function(muName.Data());
  if (!eleNorm || !muNorm || eleNorm->getVal() <= 0. || muNorm->getVal() <= 0.) {
    return false;
  }
  yield = eleNorm->getVal() + muNorm->getVal();
  source = Form("sum(%s,%s)", eleName.Data(), muName.Data());
  return true;
}

void closeFile(TFile*& file) {
  if (file) {
    file->Close();
    file = nullptr;
  }
}

}

bool SpuriousSignalTest(RooRealVar *input_mass, TString cat, TString funcType, int order, int sig, TString runPeriod, const BackgroundFitOptions& options) {
  std::cout << "[INFO] Running spurious signal test for category " << cat << " with function " << funcType << " order " << order << std::endl;
  
  // Set mass range and binning
  // double mgg_low = (cat == "VBF1" || cat == "VBF2" || cat == "VBF3") ? 105 : 100;
  double mgg_low = input_mass->getMin(), mgg_high = input_mass->getMax();
  const int bin_times = 2;
  int nBins = std::max(1, static_cast<int>(std::round((mgg_high - mgg_low) * bin_times)));
  double bin_size = nBins;
  
  // Output file path
  TString outFilePath = Form("%s/%s_%s_%s%d_%dxsig", options.outDir.c_str(),
                           runPeriod.Data(), cat.Data(), funcType.Data(), order, sig);
  ofstream output(outFilePath + ".txt", ofstream::app);
  TString inputFilePath = options.sstInputPath.c_str();
  TString inputSignalPath = options.sstSignalPath.c_str();

  TFile* fbkg = nullptr;
  TFile* fdata = nullptr;
  TFile* signalFile = nullptr;
  auto fatal = [&](const std::string& message) {
    closeFile(signalFile);
    closeFile(fdata);
    closeFile(fbkg);
    output.close();
    throw std::runtime_error(message);
  };

  TString bkgFileName = Form("%s/output_BkgMC_all.root", inputFilePath.Data());
  fbkg = TFile::Open(bkgFileName);
  if (!fbkg || fbkg->IsZombie()) {
    fatal(Form("Failed to open background MC workspace: %s", bkgFileName.Data()));
  }

  RooWorkspace* wbkg = loadWorkspace(fbkg, bkgFileName);
  TString bkgDatasetName;
  RooAbsData* bkgDataset = loadDataset(wbkg, std::vector<TString>{Form("bkgmc_13TeV_%s", cat.Data())}, bkgDatasetName);
  if (!bkgDataset) {
    fatal(Form("Background MC dataset not found for category %s", cat.Data()));
  }
  std::cout << "[INFO] Loaded background MC dataset: " << bkgDatasetName << std::endl;
  
  TString dataFileName = Form("%s/output_Data_all.root", inputFilePath.Data());
  fdata = TFile::Open(dataFileName);
  if (!fdata || fdata->IsZombie()) {
    fatal(Form("Failed to open data workspace: %s", dataFileName.Data()));
  }

  RooWorkspace* wdata = loadWorkspace(fdata, dataFileName);
  TString dataDatasetName;
  RooAbsData* dataDataset = loadDataset(
      wdata, std::vector<TString>{Form("Data_13TeV_%s", cat.Data())},
      dataDatasetName);
  if (!dataDataset) {
    fatal(Form("Data dataset not found for category %s", cat.Data()));
  }
  std::cout << "[INFO] Loaded data dataset: " << dataDatasetName << std::endl;
  
  TH1F* hbkg = datasetToHist(bkgDataset, Form("bkgmc_%s_hist", cat.Data()), nBins, mgg_low, mgg_high);
  smoothSstTemplate(hbkg);
  double mcsbevents = datasetYield(bkgDataset, mgg_low, mgg_high, options.blindLow, options.blindHigh);
  double sbevents = datasetYield(dataDataset, mgg_low, mgg_high, options.blindLow, options.blindHigh);
  
  // Scale background to match sideband data
  if (mcsbevents > 0) {
    const double bkgScale = sbevents / mcsbevents;
    std::cout << "[INFO] Sideband yields for " << cat
              << ": data = " << sbevents
              << ", bkg MC = " << mcsbevents
              << ", scale = " << bkgScale << std::endl;
    hbkg->Scale(sbevents / mcsbevents);
  } else {
    fatal(Form("Background MC sideband yield is zero for category %s", cat.Data()));
  }
  
  RooRealVar* mass = input_mass;
  mass->setRange(mgg_low, mgg_high);
  RooAbsPdf* signalPdf = nullptr;
  TString signalFileName = Form(
      "%s/outdir_combinedPDFs/CMS-HZG_combinedPDFs_%s.root",
      inputSignalPath.Data(), cat.Data());
  signalFile = TFile::Open(signalFileName);
  if (!signalFile || signalFile->IsZombie()) {
    fatal(Form("Required signal PDF file not found for category %s under %s", cat.Data(), inputSignalPath.Data()));
  }

  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(signalFile->Get("wsig_13TeV"));
  if (!ws) fatal(Form("Workspace wsig_13TeV not found in signal PDF file %s", signalFileName.Data()));

  RooRealVar* sigMass = ws->var(kMassVarName);
  RooRealVar* mh = ws->var("MH");
  signalPdf = ws->pdf(Form("combinedSigPdf_%s", cat.Data()));
  if (!sigMass || !mh || !signalPdf) {
    fatal(Form("Required signal mass variables or combined signal PDF not found in %s for category %s",
               signalFileName.Data(), cat.Data()));
  }

  mass = sigMass;
  mh->setVal(125);
  mass->setRange(mgg_low, mgg_high);

  RooArgSet* signalParams = signalPdf->getParameters(RooArgSet(*mass));
  if (signalParams) for (RooAbsArg *arg : *signalParams) {
    RooRealVar *param = dynamic_cast<RooRealVar*>(arg);
    if (!param) continue;
    if (!param->isConstant()) param->setConstant(true);
  }
  delete signalParams;
  std::cout << "[INFO] Loaded signal PDF: " << signalFileName
            << " for category " << cat << std::endl;

  std::cout << "[INFO] Spurious Signal Testing: mass range " << mass->getMin() << " - " << mass->getMax() << std::endl;

  TH1F* hsig = new TH1F(Form("sig_%s_shape", cat.Data()), Form("sig_%s_shape", cat.Data()), nBins, mgg_low, mgg_high);
  hsig->Sumw2();

  double sigevents = 0.;
  TString signalYieldSource;
  bool foundSignalYield = signalYieldFromPdfWorkspace(ws, cat, sigevents, signalYieldSource);
  if (foundSignalYield) {
    const double rawSignalYield = sigevents;
    sigevents *= kSignalPdfYieldScale;
    std::cout << "[INFO] Loaded signal yield from signal PDF file object " << signalYieldSource
              << ": raw = " << rawSignalYield
              << ", scaled = " << sigevents << " events" << std::endl;
  }

  if (!foundSignalYield || sigevents <= 0.) {
    fatal(Form("Required positive signal normalization not found in %s for category %s",
               signalFileName.Data(), cat.Data()));
  }
  std::cout << "[INFO] Total signal yield for SST: " << sigevents << std::endl;

  double bkgevents = hbkg->Integral();
  
  // Create RooDataHist objects
  RooDataHist* dbkg = new RooDataHist("bkg_mc", "dataset with x", *mass, hbkg);
  
  // Set signal and background normalization variables
  RooRealVar nsig("nsig", "nsig", sig * sigevents, -100 * sigevents, 100 * sigevents);
  RooRealVar nbkg("nbkg", "nbkg", bkgevents, 0.01 * bkgevents, 2 * bkgevents);
  
  // Create PDF model builder
  PdfModelBuilder pdfsModel;
  pdfsModel.setObsVar(mass);
  
  // Set variables to track fit results
  TString status = "Pass";
  double chi2 = 0, prob = 0, nll = 0;
  double dmc = 0, dss = 0, ss = 0, ss_mc = 0;
  double tot_err = 0, ss_cor = 0, delta = 0;
  int fit_status = 0;
  
  string typePrefix;
  // Create background model
  RooAbsPdf* bkg_model = getPdf(pdfsModel, funcType.Data(), order, &typePrefix, 
                              Form("sstest_pdf_%s_%s%d", cat.Data(), funcType.Data(), order), cat.Data(), /*doBlind*/false, options.useDiGaussian);
  if (!bkg_model) {
    std::cerr << "[ERROR] Failed to create background PDF" << std::endl;
    closeFile(signalFile);
    closeFile(fbkg);
    closeFile(fdata);
    output.close();
    return false;
  }
  
  // Fit background model to background template
  const auto weightedFitQualityPass = [](const RooFitResult* result) {
    return result && result->status() == 0 && result->covQual() >= 3 &&
           std::isfinite(result->edm()) && result->edm() <= 1.e-3;
  };
  RooFitResult* bkg_model_fit = nullptr;
  bool bkgFitPass = false;
  for (int attempt = 0; attempt < 3 && !bkgFitPass; ++attempt) {
    delete bkg_model_fit;
    bkg_model_fit = fitSstChi2Template(
        *bkg_model, *dbkg, false,
        Form("%s_%s%d_bkg_sumw2_try%d", cat.Data(), funcType.Data(), order,
             attempt + 1));
    bkgFitPass = weightedFitQualityPass(bkg_model_fit);
  }
  std::cout << "[SST_FIT] stage=background category=" << cat
            << " family=" << funcType << " order=" << order
            << " status=" << (bkg_model_fit ? bkg_model_fit->status() : -999)
            << " cov_qual=" << (bkg_model_fit ? bkg_model_fit->covQual() : -999)
            << " edm=" << (bkg_model_fit ? bkg_model_fit->edm() : -999.)
            << " passed=" << (bkgFitPass ? 1 : 0) << std::endl;
  if (!bkg_model_fit) {
    delete bkg_model;
    closeFile(signalFile);
    closeFile(fbkg);
    closeFile(fdata);
    output.close();
    return false;
  }
  if (!bkgFitPass) status = "Fail";
      
  // Plot background fit
  RooPlot* frame_bkg = mass->frame(Range(mgg_low, mgg_high), Title(Form("Background %s Fit", funcType.Data())));
  dbkg->plotOn(frame_bkg);
  bkg_model->plotOn(frame_bkg);
  bkg_model->paramOn(frame_bkg, Layout(0.34, 0.96, 0.89), Format("NEA", AutoPrecision(1)));
  bkg_model->SetName(funcType);
  
  nll = bkg_model_fit->minNll();
  int bkg_npars = bkg_model_fit->floatParsFinal().getSize();
  int bkg_ndof = bin_size - bkg_npars;
  
  output << "\t" << funcType.Data() << "\tbkg:\tnpars = " << bkg_npars 
         << " \tchi^2 = " << frame_bkg->chiSquare(bkg_npars) 
         << "\tprob = " << TMath::Prob(frame_bkg->chiSquare(bkg_npars) * bkg_ndof, bkg_ndof) 
         << "\tnll: " << nll << "\tstatus = " << bkg_model_fit->status() << std::endl;
  delete bkg_model_fit;
  
  // Create pseudo data: inject signal
  TH1F* hdata = (TH1F*)hbkg->Clone("hdata");
  hsig->Reset();
  
  // Normalize signal PDF and fill histogram
  RooAbsReal* normIntegral = signalPdf->createIntegral(*mass);
  if (!normIntegral || normIntegral->getVal() <= 0) {
    delete bkg_model;
    fatal(Form("Signal PDF normalization invalid for category %s", cat.Data()));
  }
  
  double normVal = normIntegral->getVal();
  
  // // Fill signal histogram
  for (int i = 1; i <= hsig->GetNbinsX(); i++) {
    double x = hsig->GetBinCenter(i);
    mass->setVal(x);
    double val = signalPdf->getVal() * sigevents * hsig->GetBinWidth(i) / normVal;
    if (val < 0 || val > 100) val = 0;
    hsig->SetBinContent(i, val);
  }
  
  // // Add signal to background
  hdata->Add(hsig, sig);
  RooDataHist* ddata = new RooDataHist("data_bin", "datahist with x", *mass, hdata);
  // Create combined model
  RooAddPdf* model = new RooAddPdf("model", "model", RooArgList(*signalPdf, *bkg_model), 
                                 RooArgList(nsig, nbkg));
  
  // Set canvas and plot frames
  TCanvas* canv = new TCanvas();
  RooPlot* frame_data = mass->frame(Range(mgg_low, mgg_high));
  RooPlot* frame_data_trash = mass->frame(Range(mgg_low, mgg_high));
  
  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.35);
  pad1->SetBottomMargin(0.18);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  
  // Fit combined model to data
  RooFitResult* model_fit = nullptr;
  bool sumw2FitPass = false;
  for (int attempt = 0; attempt < 3 && !sumw2FitPass; ++attempt) {
    delete model_fit;
    model_fit = fitSstChi2Template(
        *model, *ddata, true,
        Form("%s_%s%d_sb_sumw2_try%d", cat.Data(), funcType.Data(), order,
             attempt + 1));
    sumw2FitPass = weightedFitQualityPass(model_fit);
  }
  fit_status = model_fit ? model_fit->status() : -999;
  std::cout << "[SST_FIT] stage=sumw2 category=" << cat
            << " family=" << funcType << " order=" << order
            << " status=" << fit_status
            << " cov_qual=" << (model_fit ? model_fit->covQual() : -999)
            << " edm=" << (model_fit ? model_fit->edm() : -999.)
            << " passed=" << (sumw2FitPass ? 1 : 0) << std::endl;
  if (!model_fit) {
    delete bkg_model;
    closeFile(signalFile);
    closeFile(fbkg);
    closeFile(fdata);
    output.close();
    return false;
  }
  if (!sumw2FitPass) status = "Fail";
  
  // Extract signal strength and error
  ss_mc = nsig.getVal();
  dmc = nsig.getError();
  
  // Plot data and fit
  ddata->plotOn(frame_data, Name("data"), DataError(RooAbsData::SumW2));
  RooHist* plotdata = (RooHist*)frame_data->getObject(frame_data->numItems() - 1);
  
  ddata->plotOn(frame_data_trash, Name("data"), DataError(RooAbsData::SumW2));
  model->plotOn(frame_data_trash, Name("fit"));
  
  int data_npars = model_fit ? model_fit->floatParsFinal().getSize() : 1;
  int data_ndof = bin_size - data_npars;
  chi2 = frame_data_trash->chiSquare(data_npars);
  prob = TMath::Prob(chi2 * data_ndof, data_ndof);

  output << "\t" << funcType.Data() << "\tdata(MC):\tnpars = " << data_npars 
         << " \tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tnll: " << nll << "\tstatus = " << fit_status
         << std::endl;
  
  // runFit(model, ddata, &nll, &fit_status, /*max iterations*/3, /*doBlind*/false, /*cat*/cat.Data(), typePrefix, order, bin_size, model);

  // Refit using Poisson errors.
  delete model_fit;
  model_fit = nullptr;
  bool poissonFitPass = false;
  for (int attempt = 0; attempt < 3 && !poissonFitPass; ++attempt) {
    delete model_fit;
    model_fit = fitSstNllTemplate(
        *model, *ddata, true,
        Form("%s_%s%d_sb_poisson_try%d", cat.Data(), funcType.Data(), order,
             attempt + 1));
    poissonFitPass = model_fit && model_fit->status() == 0 &&
                     model_fit->covQual() >= 3 &&
                     std::isfinite(model_fit->edm()) &&
                     model_fit->edm() <= 1.e-3;
  }
  fit_status = model_fit ? model_fit->status() : -999;
  std::cout << "[SST_FIT] stage=poisson category=" << cat
            << " family=" << funcType << " order=" << order
            << " status=" << fit_status
            << " cov_qual=" << (model_fit ? model_fit->covQual() : -999)
            << " edm=" << (model_fit ? model_fit->edm() : -999.)
            << " passed=" << (poissonFitPass ? 1 : 0) << std::endl;
  if (model_fit) {
    for (int i = 0; i < model_fit->numStatusHistory(); ++i) {
      std::cout << "[SST_STATUS_HISTORY] stage=poisson category=" << cat
                << " family=" << funcType << " order=" << order
                << " index=" << i
                << " label=" << model_fit->statusLabelHistory(i)
                << " code=" << model_fit->statusCodeHistory(i) << std::endl;
    }
  }
  if (!model_fit) {
    delete bkg_model;
    closeFile(signalFile);
    closeFile(fbkg);
    closeFile(fdata);
    output.close();
    return false;
  }
  if (!poissonFitPass) status = "Fail";

  // Get signal strength and error
  ss = nsig.getVal();
  dss = nsig.getError();
  tot_err = sqrt(dss * dss + ss * ss);
  delta = fabs(ss) - 2 * dmc;
  ss_cor = (delta < 0) ? 0 : ((ss > 0) ? delta : -delta);
  
  // Test failure criteria
  const double sstThreshold = 0.2 * dss + sig * sigevents;
  if (delta > sstThreshold) status = "Fail";

  std::cout << "[SST_RESULT] category=" << cat
            << " family=" << funcType << " order=" << order
            << " injected_signal=" << sig * sigevents
            << " sumw2_signal=" << ss_mc
            << " sumw2_error=" << dmc
            << " poisson_signal=" << ss
            << " poisson_error=" << dss
            << " delta=" << delta
            << " threshold=" << sstThreshold
            << " passed=" << (status == "Pass" ? 1 : 0) << std::endl;

  output << "\t" << funcType.Data() << "\tdata(Po):\tnpars = " << data_npars 
         << " \tchi^2 = " << chi2 << "\tprob = " << prob 
         << "\tnll: " << nll << "\tstatus = " << fit_status
         << "\tsig: " << ss_mc << " +/- " << dmc 
         << "\tsig (unweighted): " << ss << " +/- " << dss 
         << "\tdelta: " << delta 
         << "\nsig_in: " << sig 
         << "\nsig_in * sigevents: " << sig * sigevents 
         << "\ttot_err: " << tot_err 
         << "\tss_cor: " << ss_cor 
         << std::endl;
  
  // Plot fit results
  model->plotOn(frame_data, Name("fit"));
  model->plotOn(frame_data, Name("background"), Components(bkg_model->GetName()), 
              LineStyle(ELineStyle::kDashed), LineColor(kGreen));
  RooCurve* nomBkgCurve = (RooCurve*)frame_data->getObject(frame_data->numItems() - 1);
  
  // Create legend
  TLegend* leg = new TLegend(0.6, 0.65, 0.88, 0.88);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->AddEntry(frame_data->findObject("data"), "Data", "ep");
  leg->AddEntry(frame_data->findObject("fit"), "Bkg+Sig", "l");
  leg->AddEntry(frame_data->findObject("background"), "Bkg", "l");
  
  // Set title and style
  model->SetName(Form("%s_model", funcType.Data()));
  frame_data->SetTitle(Form("Pseudo Data (x%d Sig), Prob: %.3f", sig, prob));
  frame_data->SetXTitle("");
  frame_data->SetLabelSize(0.042, "XY");
  frame_data->SetTitleSize(0.056, "Y");
  frame_data->SetTitleOffset(0.75, "Y");
  frame_data->Draw();
  leg->Draw("same");
  
  // Add signal curve
  signalPdf->plotOn(frame_data, Name("signal"),Normalization(ss, RooAbsReal::NumEvent), 
                  LineColor(kRed), LineWidth(4));
  
  double minResidual = 1000.;
  double maxResidual = -1000.;

  // Plot residuals and signal curve
  TGraphAsymmErrors* hdatasub = new TGraphAsymmErrors();
  int point = 0;
  for (int i = 0; i < plotdata->GetN(); i++) {
    double x, y;
    plotdata->GetPoint(i, x, y);
    double bkg = nomBkgCurve->interpolate(x);
    hdatasub->SetPoint(point, x, y - bkg);
    hdatasub->SetPointError(point, 0, 0, plotdata->GetErrorYlow(i), plotdata->GetErrorYhigh(i));
    minResidual = std::min(minResidual, y - bkg - plotdata->GetErrorYlow(i));
    maxResidual = std::max(maxResidual, y - bkg + plotdata->GetErrorYhigh(i));
    point++;
  }
  
  // Create residuals panel
  pad2->cd();
  TH1D* hdummy = new TH1D("hdummyweight", "", mgg_high - mgg_low, mgg_low, mgg_high);
  hdummy->SetStats(0);
  hdummy->SetMaximum(maxResidual+0.1*(maxResidual-minResidual));
  hdummy->SetMinimum(minResidual-0.1*(maxResidual-minResidual));
  hdummy->GetYaxis()->SetTitle("Data - Bkg PDF");
  hdummy->GetYaxis()->SetTitleOffset(0.35);
  hdummy->GetYaxis()->SetTitleSize(0.12);
  hdummy->GetYaxis()->SetLabelSize(0.09);
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  hdummy->GetXaxis()->SetTitleSize(0.12);
  hdummy->GetXaxis()->SetLabelSize(0.09);
  hdummy->Draw("HIST");
  
  hdatasub->SetMarkerStyle(8);
  hdatasub->Draw("PESAME");
  
  // Save image
  canv->SaveAs(Form("%s/%s_sstest_%s_%s%d_%dsig.png", options.outDir.c_str(),
                      cat.Data(), runPeriod.Data(), funcType.Data(), order, sig));
  if (options.plotCorr){
    plot_corr(model_fit, Form("%s/%s_sstest_corr_%s_%s%d_%dsig.png", options.outDir.c_str(),
                      cat.Data(), runPeriod.Data(), funcType.Data(), order, sig));
  }
  
  // Clean up memory
  delete model_fit;
  delete bkg_model;
  delete canv;
  delete hdummy;
  delete hdatasub;
  delete leg;
  closeFile(signalFile);
  closeFile(fbkg);
  closeFile(fdata);
  output.close();
  
  // Return test result
  return (status == "Pass");
}
