#include "../interface/WorkspaceIO.h"

#include <iostream>

#include "../interface/PlotUtils.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdfCombine.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "TKey.h"
#include "TMacro.h"
#include "TString.h"

using namespace std;

void transferMacros(TFile *inFile, TFile *outFile){
  TIter next(inFile->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    if (string(key->ReadObj()->ClassName())=="TMacro") {
      TMacro *macro = (TMacro*)inFile->Get(key->GetName());
      outFile->cd();
      macro->Write();
    }
  }
}

bool writeMultiPdfForCategory(RooWorkspace* outputws, FTestCategoryResult& result, RooRealVar* mass, RooDataSet* data, RooDataSet* dataSB, RooDataSet* dataFull, int cat, const std::string& energyExt, BackgroundFitOptions& options, bool verbose) {
  string catindexname;
  string catname;
  if (options.isFlashgg){
    catindexname = Form("pdfindex_%s_%s", options.flashggCats[cat].c_str(), energyExt.c_str());
    catname = Form("%s", options.flashggCats[cat].c_str());
  } else {
    catindexname = Form("pdfindex_%d_%s", (cat + options.catOffset), energyExt.c_str());
    catname = Form("cat%d", (cat + options.catOffset));
  }

  if (result.storedPdfs.getSize() == 0) {
    std::cerr << "[ERROR] No envelope PDFs passed GOF/KS/SST/F-test for category " << catname
              << "; skipping MultiPdf creation." << std::endl;
    return false;
  }

  RooCategory catIndex(catindexname.c_str(), "c");
  RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_hzg_%s_%s_bkgshape", catname.c_str(), energyExt.c_str()), "all pdfs", catIndex, result.storedPdfs);
  RooRealVar nBackground(Form("CMS_hzg_%s_%s_bkgshape_norm", catname.c_str(), energyExt.c_str()), "nbkg", data->sumEntries(), 0, 3 * data->sumEntries());

  std::cout << "[INFO] Check the best pdf" << std::endl;
  int bestFitPdfIndex = 0;
  if (result.storedPdfs.getSize() == 1) {
    catIndex.setIndex(0);
    std::cout << "[BEST_FIT_REUSE] category=" << catname
              << " retained_candidates=1 selected_index=0"
              << " reason=unique_candidate" << std::endl;
  } else {
    bestFitPdfIndex = getBestFitFunction(
        pdf, data, dataSB, &catIndex, options, !verbose);
  }
  catIndex.setIndex(bestFitPdfIndex, true);

  std::cout << "// ------------------------------------------------------------------------- //" << std::endl;
  std::cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs" << std::endl;
  result.storedPdfs.Print();
  std::cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << result.storedPdfs.at(bestFitPdfIndex)->GetName() << std::endl;
  std::cout << "// ------------------------------------------------------------------------- //" << std::endl;
  std::cout << "[INFO] Simple check of index " << result.simpleBestFitPdfIndex << std::endl;

  mass->setBins(options.nBinsForMass);
  RooDataHist dataBinned(Form("roohist_data_mass_%s", catname.c_str()), "data", *mass, *dataFull);

  std::cout << "[INFO] Import pdf " << pdf->GetName() << std::endl;
  outputws->import(*pdf);
  std::cout << "[INFO] Import nBackground " << nBackground.GetName() << std::endl;
  outputws->import(nBackground);
  std::cout << "[INFO] Import dataBinned " << dataBinned.GetName() << std::endl;
  outputws->import(dataBinned);
  plot(mass, pdf, &catIndex, data, Form("%s/multipdf_%s", options.outDir.c_str(), catname.c_str()), options.flashggCats, cat, bestFitPdfIndex, options);
  return true;
}
