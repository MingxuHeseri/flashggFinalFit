#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/program_options.hpp"

#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooNumIntConfig.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TError.h"
#include "TFile.h"

#include <filesystem>
#include "../interface/BackgroundFitOptions.h"
#include "../interface/FTest.h"
#include "../interface/PdfFactory.h"
#include "../interface/PlotUtils.h"
#include "../interface/WorkspaceIO.h"

using namespace std;
namespace po = boost::program_options;

// Persist the requested accuracy, including on dependency functions, so a
// fresh downstream process does not silently return to ROOT's defaults.
static void persistIntegrationPrecision(RooWorkspace& workspace, double epsilon) {
  if (epsilon == 0.) return;
  RooArgSet nodes(workspace.allPdfs());
  nodes.add(workspace.allFunctions(), true);
  for (RooAbsArg* node : nodes) {
    auto* real = dynamic_cast<RooAbsReal*>(node);
    if (!real) continue;
    RooNumIntConfig config(*real->getIntegratorConfig());
    config.setEpsAbs(epsilon);
    config.setEpsRel(epsilon);
    real->setIntegratorConfig(config);
  }
}

int main(int argc, char* argv[]){
  configureBackgroundPlotStyle();

  BackgroundFitOptions options;
  double integrationEpsilon = 0.;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", "Show help")
    ("infilename,i", po::value<string>(&options.fileName), "In file name")
    ("ncats,c", po::value<int>(&options.nCats)->default_value(5), "Number of categories")
    ("singleCat", po::value<int>(&options.singleCategory)->default_value(-1), "Run A single Category")
    ("datfile,d", po::value<string>(&options.datfile)->default_value("dat/fTest.dat"), "Right results to datfile for BiasStudy")
    ("outDir,D", po::value<string>(&options.outDir)->default_value("plots/fTest"), "Out directory for plots")
    ("saveMultiPdf", po::value<string>(&options.outfilename), "Save a MultiPdf model with the appropriate pdfs")
    ("funclist", po::value<string>(&options.functionList)->default_value(options.functionList), "Function list for bkg fitting")
    ("useDiGaussian", "Use the di-Gaussian convolution for step x Gaussian families")
    ("runFtestCheckWithToys", "When running the F-test, use toys to calculate pvals (and make plots) ")
    ("is2011", "Run 2011 config")
    ("is2012", "Run 2012 config")
    ("unblind", "Do not blind plots")
    ("blindFit", "Blind fits in signal region (120-130 GeV)")
    ("plot_corr", "Plot correlation matrix of bkg functions")
    ("nsig_in", po::value<int>(&options.nsigIn)->default_value(0), "Inject signal in spurious signal test")
    ("skipSSTest", "Skip SpuriousSignalTest in the candidate selection chain")
    ("objectiveOnly", "Run the requested candidate fits and objective report without creating a MultiPdf")
    ("skipCrossObjectiveDiagnostics", po::bool_switch(&options.skipCrossObjectiveDiagnostics), "Skip redundant post-fit sideband NLL/chi2 reporting; preserve fit and all candidate acceptance checks")
    ("sstInputPath", po::value<string>(&options.sstInputPath)->default_value(options.sstInputPath), "Directory containing SST background MC/data workspaces")
    ("sstSignalPath", po::value<string>(&options.sstSignalPath)->default_value(options.sstSignalPath), "Directory or ROOT file containing SST signal combined PDFs")
    ("dochi2Fit", "Use Chi2 fit instead of NLL fit")
    ("minPdfOrder", po::value<int>(&options.minPdfOrder)->default_value(options.minPdfOrder), "Minimum PDF order scanned by the envelope F-test")
    ("maxPdfOrder", po::value<int>(&options.maxPdfOrder)->default_value(options.maxPdfOrder), "Maximum PDF order scanned by the envelope F-test")
    ("candidateOrders", po::value<string>(&options.candidateOrdersStr)->default_value(""), "Exact comma-separated PDF orders to scan; empty uses the inclusive min/max range")
    ("isFlashgg", po::value<int>(&options.isFlashgg)->default_value(1), "Use Flashgg output ")
    ("isData", po::value<bool>(&options.isData)->default_value(0), "Use Data not MC ")
    ("flashggCats,f", po::value<string>(&options.flashggCatsStr)->default_value("UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,TTHHadronicTag,TTHLeptonicTag,VHHadronicTag,VHTightTag,VHLooseTag,VHEtTag"), "Flashgg category names to consider")
    ("year", po::value<string>(&options.year)->default_value("2016"), "Dataset year")
    ("catOffset", po::value<int>(&options.catOffset)->default_value(0), "Category numbering scheme offset")
    ("mgg_low", po::value<double>(&options.mggLow)->default_value(95), "Lower bound for mgg")
    ("mgg_high", po::value<double>(&options.mggHigh)->default_value(170), "Upper bound for mgg")
    ("minimizerEps", po::value<double>(&options.minimizerEps)->default_value(options.minimizerEps), "RooMinimizer epsilon for explicit blinded sideband fits")
    ("integrationEps", po::value<double>(&integrationEpsilon)->default_value(0.), "Absolute and relative numerical integration accuracy, persisted in output; zero preserves ROOT defaults")
    ("minimizerStrategy", po::value<int>(&options.minimizerStrategy)->default_value(options.minimizerStrategy), "Minuit strategy for explicit blinded sideband fits (0, 1, or 2)")
    ("minimumGofProbability", po::value<double>(&options.minimumGofProbability)->default_value(options.minimumGofProbability), "Minimum exclusive Pearson and Poisson sideband GOF probability")
    ("minimumKsProbability", po::value<double>(&options.minimumKsProbability)->default_value(options.minimumKsProbability), "Minimum exclusive sideband KS probability")
    ("minimumCovQual", po::value<int>(&options.minimumCovQual)->default_value(options.minimumCovQual), "Minimum covariance quality required for candidate selection")
    ("maximumEdm", po::value<double>(&options.maximumEdm)->default_value(options.maximumEdm), "Maximum EDM required for candidate selection")
    ("requireFitStatusZero", "Require Minuit fit status zero for candidate selection")
    ("skipFitQualityGate", "Record Minuit status, covariance quality, EDM, and boundary diagnostics without using them as candidate-selection vetoes")
    ("gofBinsPerGeV", po::value<double>(&options.gofBinsPerGeV)->default_value(options.gofBinsPerGeV), "Bins per GeV for sideband GOF and KS validation")
    ("sidebandChi2Mode", po::value<string>(&options.sidebandChi2Mode)->default_value(options.sidebandChi2Mode), "Blinded bin expectations: integral (default) or midpoint")
    ("gofIntegralMode", po::value<string>(&options.gofIntegralMode)->default_value(options.gofIntegralMode), "GOF bin integration: auto (shared raw normalization for ExpBernstein/ExpChebychev erf PDFs), raw (explicit), or normalized (legacy)")
    ("plotNormalizationMode", po::value<string>(&options.plotNormalizationMode)->default_value(options.plotNormalizationMode), "Plot normalization: auto (direct sideband fraction for validated erf PDFs), or legacy (RooFit projection integral)")
    ("verbose,v", "Run with more output")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  options.updateBinningFromMassRange();
  if (options.sidebandChi2Mode != "midpoint" && options.sidebandChi2Mode != "integral") {
    std::cerr << "[ERROR] --sidebandChi2Mode must be integral or midpoint" << std::endl;
    return 1;
  }
  if (options.plotNormalizationMode != "auto" && options.plotNormalizationMode != "legacy") {
    std::cerr << "[ERROR] --plotNormalizationMode must be auto or legacy" << std::endl;
    return 1;
  }
  if (options.gofIntegralMode != "auto" && options.gofIntegralMode != "raw" &&
      options.gofIntegralMode != "normalized") {
    std::cerr << "[ERROR] --gofIntegralMode must be auto, raw, or normalized" << std::endl;
    return 1;
  }

  if (!std::isfinite(integrationEpsilon) || integrationEpsilon < 0. || integrationEpsilon >= 1.) {
    std::cerr << "[ERROR] --integrationEps must be zero or finite in (0, 1)" << std::endl;
    return 1;
  }
  if (integrationEpsilon > 0.) {
    RooAbsReal::defaultIntegratorConfig()->setEpsAbs(integrationEpsilon);
    RooAbsReal::defaultIntegratorConfig()->setEpsRel(integrationEpsilon);
    std::cout << "[INTEGRATION_CONFIG] eps_abs=" << integrationEpsilon
              << " eps_rel=" << integrationEpsilon << " persist=1" << std::endl;
  }

  if (options.minimizerEps <= 0.) {
    std::cerr << "[ERROR] --minimizerEps must be positive" << std::endl;
    return 1;
  }
  if (options.minimizerStrategy < 0 || options.minimizerStrategy > 2) {
    std::cerr << "[ERROR] --minimizerStrategy must be 0, 1, or 2" << std::endl;
    return 1;
  }
  if (options.gofBinsPerGeV <= 0.) {
    std::cerr << "[ERROR] --gofBinsPerGeV must be positive" << std::endl;
    return 1;
  }
  if (options.minimumGofProbability < 0. || options.minimumGofProbability >= 1.) {
    std::cerr << "[ERROR] --minimumGofProbability must be in [0, 1)" << std::endl;
    return 1;
  }
  if (options.minimumKsProbability < 0. || options.minimumKsProbability >= 1.) {
    std::cerr << "[ERROR] --minimumKsProbability must be in [0, 1)" << std::endl;
    return 1;
  }
  if (options.maxPdfOrder < 0 || options.maxPdfOrder > 20) {
    std::cerr << "[ERROR] --maxPdfOrder must be between 0 and 20" << std::endl;
    return 1;
  }
  if (options.minPdfOrder < 0 || options.minPdfOrder > options.maxPdfOrder) {
    std::cerr << "[ERROR] --minPdfOrder must be between 0 and --maxPdfOrder" << std::endl;
    return 1;
  }
  if (options.minimumCovQual < 0 || options.minimumCovQual > 3) {
    std::cerr << "[ERROR] --minimumCovQual must be between 0 and 3" << std::endl;
    return 1;
  }
  if (!(options.maximumEdm > 0.)) {
    std::cerr << "[ERROR] --maximumEdm must be positive" << std::endl;
    return 1;
  }
  if (!options.candidateOrdersStr.empty()) {
    std::vector<std::string> orderTokens;
    boost::split(orderTokens, options.candidateOrdersStr, boost::is_any_of(","));
    for (const std::string& token : orderTokens) {
      try {
        const int candidateOrder = std::stoi(token);
        if (candidateOrder < options.minPdfOrder || candidateOrder > options.maxPdfOrder) {
          throw std::out_of_range("candidate order outside min/max range");
        }
        options.candidateOrders.push_back(candidateOrder);
      } catch (const std::exception&) {
        std::cerr << "[ERROR] Invalid --candidateOrders entry: " << token << std::endl;
        return 1;
      }
    }
    std::sort(options.candidateOrders.begin(), options.candidateOrders.end());
    options.candidateOrders.erase(std::unique(options.candidateOrders.begin(), options.candidateOrders.end()), options.candidateOrders.end());
  }

  if (vm.count("help")) { cout << desc << endl; return 0; }
  if (vm.count("is2011")) options.is2011 = true;
  if (vm.count("unblind")) options.blind = false;
  if (vm.count("blindFit")) options.blindFit = true;
  if (vm.count("dochi2Fit")) options.doChi2Fit = true;
  if (vm.count("plot_corr")) options.plotCorr = true;
  if (vm.count("useDiGaussian")) options.useDiGaussian = true;
  if (vm.count("verbose")) options.verbose = true;
  if (vm.count("runFtestCheckWithToys")) options.runFtestCheckWithToys = true;
  if (vm.count("skipSSTest")) options.skipSSTest = true;
  if (vm.count("objectiveOnly")) options.objectiveOnly = true;
  if (vm.count("requireFitStatusZero")) options.requireFitStatusZero = true;
  if (vm.count("skipFitQualityGate")) options.gateFitQuality = false;
  options.saveMultiPdf = vm.count("saveMultiPdf");
  if (options.useDiGaussian) std::cout << "[INFO] Di-Gaussian option enabled for step x Gaussian families" << std::endl;

  if (!options.skipSSTest && (options.saveMultiPdf || options.objectiveOnly) &&
      (options.sstInputPath.empty() || options.sstSignalPath.empty())) {
    std::cerr << "[ERROR] Supply --sstInputPath and --sstSignalPath for the template test, or use --skipSSTest" << std::endl;
    return 1;
  }
  FILE *dfile = fopen(options.datfile.c_str(), "a");
  if (!dfile) {
    std::cerr << "[ERROR] Failed to open dat file for writing: " << options.datfile << std::endl;
    return 1;
  }

  if (!options.verbose) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    RooMsgService::instance().setSilentMode(true);
    gErrorIgnoreLevel = kWarning;
  }

  boost::split(options.flashggCats, options.flashggCatsStr, boost::is_any_of(","));
  std::cout << "[INFO] SaveMultiPdf? " << options.saveMultiPdf << std::endl;

  int startingCategory = 0;
  if (options.isFlashgg == 1) options.nCats = options.flashggCats.size();
  if (options.singleCategory > -1) {
    if (options.singleCategory >= options.nCats) {
      std::cerr << "[ERROR] Requested singleCat " << options.singleCategory
                << " but only " << options.nCats << " categories are configured" << std::endl;
      return 1;
    }
    options.nCats = options.singleCategory + 1;
    startingCategory = options.singleCategory;
  }

  if (options.verbose) std::cout << "[INFO] SaveMultiPdf? " << options.saveMultiPdf << std::endl;
  if (options.verbose) std::cout << "[INFO] Sideband Fit? " << options.blindFit << std::endl;

  TFile *outputfile = nullptr;
  RooWorkspace *outputws = nullptr;
  if (options.saveMultiPdf) {
    outputfile = new TFile(options.outfilename.c_str(), "RECREATE");
    outputws = new RooWorkspace();
    outputws->SetName("multipdf");
  }

  std::filesystem::create_directories(options.outDir);
  TFile *inFile = TFile::Open(options.fileName.c_str());
  RooWorkspace *inWS = inFile ? (RooWorkspace*)inFile->Get("CMS_hzg_workspace") : nullptr;
  if (options.verbose) std::cout << "[INFO]  inWS open " << inWS << std::endl;
  if (!inWS) {
    std::cerr << "[ERROR] Cannot find workspace CMS_hzg_workspace in " << options.fileName << std::endl;
    return 1;
  }

  if (options.saveMultiPdf) {
    transferMacros(inFile, outputfile);
    outputws->import(backgroundIntLumi());
    outputws->import(backgroundSqrtS());
    std::cout << "[INFO] got intL and sqrts " << &backgroundIntLumi() << ", " << &backgroundSqrtS() << std::endl;
  }

  vector<string> functionClasses;
  boost::split(functionClasses, options.functionList, boost::is_any_of(","));
  std::cout << "[INFO] Function list: ";
  for (const auto& s : functionClasses) std::cout << s << " ";
  std::cout << "\n";

  map<string, string> namingMap = defaultNamingMap();

  FILE *resFile = nullptr;

  vector<map<string, int> > choicesVec;
  vector<map<string, vector<int> > > choicesEnvelopeVec;

  PdfModelBuilder pdfsModel;
  RooRealVar *mass = (RooRealVar*)inWS->var("CMS_hzg_mass");
  mass->setRange(options.mggLow, options.mggHigh);
  mass->setBins(options.nBinsForMass);
  std::cout << "[INFO] Got mass from ws " << mass << std::endl;
  pdfsModel.setObsVar(mass);

  std::string energyExt = options.is2011 ? "7TeV" : "8TeV";
  if (options.isFlashgg) {
    if (options.year.find("run2") != std::string::npos) energyExt = "13TeV";
    else if (options.year.find("run3") != std::string::npos) energyExt = "13.6TeV";
    else if (options.year.find("all") != std::string::npos) energyExt = "13_13p6TeV";
    else energyExt = Form("%s_13TeV", options.year.c_str());
  }

  bool multiPdfWriteFailed = false;

  for (int cat = startingCategory; cat < options.nCats; cat++) {
    string catname = options.isFlashgg ? Form("%s", options.flashggCats[cat].c_str()) : Form("cat%d", cat);
    fprintf(dfile, "\ncat=%d\n", (cat + options.catOffset));

    RooDataSet *dataFull0 = nullptr;
    if (options.isData) {
      dataFull0 = (RooDataSet*)inWS->data(Form("Data_13TeV_%s", catname.c_str()));
      if (options.verbose) std::cout << "[INFO] opened data for  " << Form("Data_%s", catname.c_str()) << " - " << dataFull0 << std::endl;
    } else {
      dataFull0 = (RooDataSet*)inWS->data(Form("data_mass_%s", catname.c_str()));
      if (options.verbose) std::cout << "[INFO] opened data for  " << Form("data_mass_%s", catname.c_str()) << " - " << dataFull0 << std::endl;
    }
    if (!dataFull0) {
      std::cerr << "[ERROR] Cannot find input dataset for category " << catname << std::endl;
      return 1;
    }

    RooDataSet *dataFull = (RooDataSet*)dataFull0->emptyClone();
    for (int i = 0; i < dataFull0->numEntries(); i++) {
      double m = dataFull0->get(i)->getRealValue("CMS_hzg_mass");
      if ((m < (options.mggLow + 0.01)) or (m > (options.mggHigh - 0.01))) continue;
      dataFull->add(*dataFull0->get(), 1.0);
    }

    RooDataSet *dataSB = (RooDataSet*)dataFull0->emptyClone();
    for (int i = 0; i < dataFull0->numEntries(); i++) {
      double m = dataFull0->get(i)->getRealValue("CMS_hzg_mass");
      if ((m < (options.mggLow + 0.01)) or (m > (options.mggHigh - 0.01))) continue;
      if ((m < options.blindLow) or (m > options.blindHigh)) dataSB->add(*dataFull0->get(), 1.0);
    }

    string thisdataBinnedName;
    string thisdataSBBinnedName;
    if (options.isFlashgg) {
      thisdataBinnedName = Form("roohist_data_mass_%s", options.flashggCats[cat].c_str());
      thisdataSBBinnedName = Form("roohist_dataSB_mass_%s", options.flashggCats[cat].c_str());
    } else {
      thisdataBinnedName = Form("roohist_data_mass_cat%d", cat);
      thisdataSBBinnedName = Form("roohist_dataSB_mass_cat%d", cat);
    }

    // Keep the event-level datasets for the likelihood fit.  The binned
    // RooDataHist is constructed explicitly only where a binned observable
    // (chi2 fit, GOF, or workspace output) is required.  Casting a
    // RooDataHist to RooDataSet here made the nominal NLL fit binned and was
    // undefined behaviour.
    RooDataSet *data = dataFull;
    RooDataSet *dataSBForFit = dataSB;

    TString runPeriod = "all";
    if (options.year.find("run3") != std::string::npos) runPeriod = "run3";
    else if (options.year.find("run2") != std::string::npos) runPeriod = "run2";

    FTestCategoryInputs ftestInput;
    ftestInput.cat = cat;
    ftestInput.catOffset = options.catOffset;
    ftestInput.catName = catname;
    ftestInput.runPeriod = runPeriod;
    ftestInput.energyExt = energyExt;
    ftestInput.functionClasses = functionClasses;
    ftestInput.mass = mass;
    ftestInput.data = data;
    ftestInput.dataSB = dataSBForFit;
    ftestInput.dataFull = dataFull;
    ftestInput.pdfsModel = &pdfsModel;
    ftestInput.dfile = dfile;
    ftestInput.resFile = resFile;

    FTestCategoryResult ftestResult = runEnvelopeFTest(ftestInput, options);

    choicesVec.push_back(ftestResult.choices);
    choicesEnvelopeVec.push_back(ftestResult.choicesEnvelope);

    if (options.saveMultiPdf) {
      if (!writeMultiPdfForCategory(outputws, ftestResult, mass, data, dataSBForFit, dataFull, cat, energyExt, options, options.verbose)) {
        multiPdfWriteFailed = true;
        break;
      }
    }
  }

  if (multiPdfWriteFailed) {
    if (options.saveMultiPdf) {
      persistIntegrationPrecision(*outputws, integrationEpsilon);
      outputfile->cd();
      outputws->Write();
      outputfile->Close();
    }
    inFile->Close();
    if (dfile) fclose(dfile);
    std::cerr << "[ERROR] No valid envelope PDFs for at least one category; no complete MultiPdf was produced." << std::endl;
    return 2;
  }

  if (options.saveMultiPdf) {
    persistIntegrationPrecision(*outputws, integrationEpsilon);
    outputfile->cd();
    outputws->Write();
    outputfile->Close();
  }

  cout << "[RESULT] Recommended options" << endl;
  for (int cat = startingCategory; cat < options.nCats; cat++) {
    cout << "Cat " << cat << endl;
    fprintf(dfile, "\ncat=%d (code order, not true order)\n", (cat + options.catOffset));
    for (map<string, int>::iterator it = choicesVec[cat-startingCategory].begin(); it != choicesVec[cat-startingCategory].end(); it++) {
      cout << "\t" << it->first << " - " << it->second << endl;
      fprintf(dfile, "choices_vec(truth)=%s:%d:%s%d\n", it->first.c_str(), it->second, namingMap[it->first].c_str(), it->second);
    }
    for (map<string, vector<int> >::iterator it = choicesEnvelopeVec[cat-startingCategory].begin(); it != choicesEnvelopeVec[cat-startingCategory].end(); it++) {
      vector<int> ords = it->second;
      for (vector<int>::iterator ordit = ords.begin(); ordit != ords.end(); ordit++) {
        fprintf(dfile, "choices_envelope_vec=%s:%d:%s%d\n", it->first.c_str(), *ordit, namingMap[it->first].c_str(), *ordit);
      }
    }
    fprintf(dfile, "\n");
  }

  inFile->Close();
  if (dfile) fclose(dfile);
  cout << "[INFO] Program completed successfully" << endl;
  return 0;
}
