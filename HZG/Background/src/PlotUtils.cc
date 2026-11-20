#include "../interface/PlotUtils.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>

#include "../interface/GoodnessOfFit.h"
#include "../interface/PlotNormalization.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "TBox.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH2.h"
#include "TIterator.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"

#include "../interface/PlotStyle.h"


using namespace std;
using namespace RooFit;

namespace {
void plotBackgroundCurve(RooRealVar& mass, RooAbsPdf& pdf, RooPlot& frame,
                         const BackgroundFitOptions& options,
                         int color = kBlue, int style = 1) {
  const auto start = std::chrono::steady_clock::now();
  const bool fast = options.plotNormalizationMode == "auto" &&
      plotWithSharedNormalization(pdf, mass, frame, color, style);
  if (!fast) pdf.plotOn(&frame, LineColor(color), LineStyle(style));
  std::cout << "[PLOT_NORMALIZATION] pdf=" << pdf.GetName()
            << " requested=" << options.plotNormalizationMode
            << " backend=" << (fast ? "direct_sideband_fraction" : "legacy")
            << " seconds=" << std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - start).count() << std::endl;
}
}

void configureBackgroundPlotStyle() {
  setTDRStyle();
  writeExtraText = true;
  extraText = "Preliminary";
  lumi_8TeV = "19.1 fb^{-1}";
  lumi_7TeV = "4.9 fb^{-1}";
  lumi_sqrtS = "13 TeV + 13.6 TeV";
}

void plot(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name,vector<string> flashggCats_, int status, double *prob, BackgroundFitOptions& options){
  
  std::cout << "[INFO] Start plotting " << name.c_str() << std::endl;
  std::cout << "[INFO] Plotting mass range: " << mass->getMin() << " - " << mass->getMax() << std::endl;

  const int mgg_low = mass->getMin();
  const int mgg_high = mass->getMax();
  const int nBinsForMass = std::max(
      1, static_cast<int>(std::round(
             (mgg_high - mgg_low) * options.gofBinsPerGeV)));
  double gofChi2 = 0.;
  int gofNdof = 0;
  *prob = getGoodnessOfFit(mass,pdf,data,name, options, options.blindFit, &gofChi2, &gofNdof);
  RooPlot *plot = mass->frame();
  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (options.blind) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(nBinsForMass));

 // data->plotOn(plot,Binning(mgg_high-mgg_low));
  TCanvas *canv = new TCanvas();
  plotBackgroundCurve(*mass, *pdf, *plot, options);
  pdf->paramOn(plot, RooFit::Layout(0.34, 0.96, 0.89), RooFit::Format("NEA", AutoPrecision(1)));
  plot->getAttText()->SetTextSize(0.025);
  if (options.blind) plot->SetMinimum(0.0001);
  plot->SetTitle("");
  plot->Draw();
  TLatex *lat = new TLatex();
  lat->SetNDC();
  lat->SetTextFont(42);
  lat->SetTextSize(0.028);
  lat->DrawLatex(0.14,0.84,Form("%sfit quality: #chi^{2}=%.3f, ndf=%d, p=%.2f, status=%d", options.blindFit ? "Sideband " : "", gofChi2, gofNdof, *prob, status));
  CMS_lumi( canv, 0, 0, "");
  canv->SaveAs(name.c_str());
 	
	//plot_chi2->Draw();
  //canv->SaveAs((name+"debug").c_str());

  delete plot;
  delete canv;
  delete lat;
}

void plot(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf, const BackgroundFitOptions& options){
  
  std::cout << "[INFO] Start plotting " << name.c_str() << std::endl;
  const int mgg_low = mass->getMin();
  const int mgg_high = mass->getMax();
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TLegend *leg = new TLegend(0.5,0.55,0.92,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (options.blind) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(4*(mgg_high-mgg_low))); 
  TCanvas *canv = new TCanvas();
  ///start extra bit for ratio plot///
  RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
  TPad *pad1 = new TPad("pad1","pad1",0,0.25,1,1);
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.35);
  pad1->SetBottomMargin(0.18);
  pad2->SetTopMargin(0.00001);
  pad2->SetBottomMargin(0.25);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  // enf extra bit for ratio plot///

  int currentIndex = catIndex->getIndex();
  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
  leg->AddEntry(datLeg,Form("Data sidebands - %s",flashggCats_[cat].c_str()),"LEP");
  int style=1;
  RooAbsPdf *pdf;
  RooCurve *nomBkgCurve;
  int bestcol= -1;
  for (int icat=0;icat<catIndex->numTypes();icat++){
    int col;
    if (icat<=6) col=color[icat];
    else {col=kBlack; style++;}
    catIndex->setIndex(icat);
    // pdfs->getCurrentPdf()->fitTo(*data,RooFit::Minos(0),RooFit::Minimizer("Minuit2","minimize"),RooFit::SumW2Error(kTRUE));	 //FIXME
    plotBackgroundCurve(*mass, *pdfs->getCurrentPdf(), *plot, options, col, style);
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    if (bestFitPdf==icat) {
    ext=" (best fit)";
    pdf= pdfs->getCurrentPdf();
    nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);
    bestcol = col;
    }
    leg->AddEntry(pdfLeg,Form("%s%s",pdfs->getCurrentPdf()->GetName(),ext.c_str()),"L");
  }
  plot->SetTitle(Form("Category %s",flashggCats_[cat].c_str()));
  if (options.blind) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw("same");
  CMS_lumi( pad1, 0, 0, "");
  ///start extra bit for ratio plot///
  TH1D *hbplottmp = (TH1D*) pdf->createHistogram("hbplottmp",*mass,Binning(mgg_high-mgg_low,mgg_low,mgg_high));
  hbplottmp->Scale(plotdata->Integral());
  hbplottmp->Draw("same");
  int npoints = plotdata->GetN();
  double xtmp,ytmp;//
  int point =0;
  TGraphAsymmErrors *hdataratio = new TGraphAsymmErrors(npoints);
  //hdatasub->SetMarkerSize(defmarkersize);
  for (int ipoint=0; ipoint<npoints; ++ipoint) {
  //double bkgval = hbplottmp->GetBinContent(ipoint+1);
  plotdata->GetPoint(ipoint, xtmp,ytmp);
  double bkgval = nomBkgCurve->interpolate(xtmp);
  if (options.blind) {
   if ((xtmp > 120 ) && ( xtmp < 130) ) continue;
  }
  // std::cout << "[INFO] plotdata->Integral() " <<  plotdata->Integral() << " ( bins " << npoints  << ") hbkgplots[i]->Integral() " << hbplottmp->Integral() << " (bins " << hbplottmp->GetNbinsX() << std::endl;
 if (bkgval <= 0.) continue;
 double errhi = plotdata->GetErrorYhigh(ipoint) / bkgval;
 double errlow = plotdata->GetErrorYlow(ipoint) / bkgval;
       
 //std::cout << "[INFO]  Channel " << name  << " errhi " << errhi << " errlow " << errlow  << std::endl;
 //std::cout << "[INFO] Channel  " << name << " setting point " << point <<" : xtmp "<< xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << std::endl;
 bool drawZeroBins_ =1;
 if (!drawZeroBins_) if(fabs(ytmp)<1e-5) continue; 
 hdataratio->SetPoint(point,xtmp,ytmp/bkgval);
 hdataratio->SetPointError(point,0.,0.,errlow,errhi );
 point++;
  } 
  pad2->cd();
  TH1 *hdummy = new TH1D("hdummyweight","",mgg_high-mgg_low,mgg_low,mgg_high);
  hdummy->SetMaximum(point > 0 ? std::min(2.0, std::max(1.2, hdataratio->GetHistogram()->GetMaximum()+0.15)) : 1.2);
  hdummy->SetMinimum(point > 0 ? std::max(0.0, std::min(0.8, hdataratio->GetHistogram()->GetMinimum()-0.15)) : 0.8);
  hdummy->GetYaxis()->SetTitle("Data / best-fit");
  hdummy->GetYaxis()->SetTitleSize(0.12);
  hdummy->GetXaxis()->SetTitle("m_{ll#gamma} (GeV)");
  hdummy->GetXaxis()->SetTitleSize(0.12);
  hdummy->Draw("HIST");
  hdummy->GetYaxis()->SetNdivisions(808);

  TLine *line3 = new TLine(mgg_low,1.,mgg_high,1.);
  line3->SetLineColor(bestcol);
  //line3->SetLineStyle(kDashed);
  line3->SetLineWidth(5.0);
  line3->Draw();
  hdataratio->Draw("PESAME");
  // enf extra bit for ratio plot///
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  catIndex->setIndex(currentIndex);
  delete canv;
}

void plot(RooRealVar *mass, map<string,RooAbsPdf*> pdfs, RooDataSet *data, string name, vector<string> flashggCats_, int cat, int bestFitPdf, const BackgroundFitOptions& options){
  
  std::cout << "[INFO] Start plotting " << name.c_str() << std::endl;
  const int mgg_low = mass->getMin();
  const int mgg_high = mass->getMax();
  int color[7] = {kBlue,kRed,kMagenta,kGreen+1,kOrange+7,kAzure+10,kBlack};
  TCanvas *canv = new TCanvas();
  TLegend *leg = new TLegend(0.6,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  RooPlot *plot = mass->frame();

  mass->setRange("unblindReg_1",mgg_low,120);
  mass->setRange("unblindReg_2",130,mgg_high);
  if (options.blind) {
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_1"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),CutRange("unblindReg_2"));
    data->plotOn(plot,Binning(mgg_high-mgg_low),Invisible());
  }
  else data->plotOn(plot,Binning(4*(mgg_high-mgg_low)));

  TObject *datLeg = plot->getObject(int(plot->numItems()-1));
	if(flashggCats_.size() >0){
  leg->AddEntry(datLeg,Form("Data sidebands - %s",flashggCats_[cat].c_str()),"LEP");
	} else {
  leg->AddEntry(datLeg,Form("Data sidebands - %d",cat),"LEP");
	}
  int i=0;
  int style=1;
  for (map<string,RooAbsPdf*>::iterator it=pdfs.begin(); it!=pdfs.end(); it++){
    int col;
    if (i<=6) col=color[i];
    else {col=kBlack; style++;}
    cout << "[INFO] Plotting " << it->first << " with color " << col << endl;
    plotBackgroundCurve(*mass, *it->second, *plot, options, col, style);
    TObject *pdfLeg = plot->getObject(int(plot->numItems()-1));
    std::string ext = "";
    // if (bestFitPdf==i) ext=" (Best Fit Pdf) ";
    leg->AddEntry(pdfLeg,Form("Fit model: %s%s",it->first.c_str(),ext.c_str()),"L");
    i++;
  }
  plot->SetTitle(Form(" %s",flashggCats_[cat].c_str()));
  if (options.blind) plot->SetMinimum(0.0001);
  plot->Draw();
  leg->Draw();
  CMS_lumi( canv, 0, 0, "");
  canv->SaveAs(Form("%s.pdf",name.c_str()));
  canv->SaveAs(Form("%s.png",name.c_str()));
  delete canv;
}

void plot_corr(RooFitResult* model_fit, TString outputPath){
  
  TCanvas* c_corr = new TCanvas("c_corr", "Correlations", 800, 600);
  TH2* hCorr = model_fit->correlationHist();
  hCorr->SetMinimum(-1.0);
  hCorr->SetMaximum(1.0);
  hCorr->SetMarkerSize(1.5);
  hCorr->SetMarkerColor(kWhite);
  hCorr->Draw("COLZ TEXT");
  c_corr->Update();
  RooArgList params_corr = model_fit->floatParsFinal();
  TString shortNames[params_corr.getSize()];
  for (int i = 0; i < params_corr.getSize(); i++) {
    TString fullName = params_corr[i].GetName();
    TObjArray* parts = fullName.Tokenize("_");
    if (parts->GetEntries() > 0) {
      // 获取最后一个部分
      shortNames[i] = ((TObjString*)parts->At(parts->GetEntries()-1))->GetString();
    } else {
      shortNames[i] = fullName;
    }
    
    // 设置轴标签
    hCorr->GetXaxis()->SetBinLabel(i+1, shortNames[i]);
    hCorr->GetYaxis()->SetBinLabel(params_corr.getSize() - i, shortNames[i]);
  }
  c_corr->SaveAs(outputPath);
  delete c_corr;

}
