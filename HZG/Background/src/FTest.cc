#include "../interface/FTest.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>

#include "../interface/FitUtils.h"
#include "../interface/GoodnessOfFit.h"
#include "../interface/PdfFactory.h"
#include "../interface/PlotUtils.h"
#include "../interface/SpuriousSignalTest.h"
#include "TMath.h"

using namespace std;

FTestCategoryResult::FTestCategoryResult() : storedPdfs("store") {}

int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooDataSet *data_SB, RooCategory *cat, const BackgroundFitOptions& options, bool silent){
	double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indices = cat->numTypes();
	string typePrefix = ""; // 添加typePrefix变量声明
	int order = 0; // 添加order变量声明
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		//params->Print("V");
	}
	
	//bkg->setDirtyInhibit(1);
	//RooAbsReal *nllm = bkg->createNLL(*data);
	//RooMinimizer minim(*nllm);
	//minim.setStrategy(1);
	
	for (int id=0;id<number_of_indices;id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		// if (!silent) {
			/*
			std::cout << "BEFORE  MAKING FIT" << std::endl;
			params->Print("V");
			std::cout << "-----------------------" << std::endl;		
			*/
		// }
		
		//minim.minimize("Minuit2","minimize");
		double minNll=0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus=1;		
		runFit(bkg->getCurrentPdf(), data, data_SB, &minNll, &fitStatus, options, options.blindFit,
      /*cat*/cat->GetName(), typePrefix, order, options.nBinsForMass, bkg->getCurrentPdf(), options.doChi2Fit);

		// Add the penalty
		minNll=minNll+bkg->getCorrection();

		if (!silent) {
			std::cout << "[INFO] AFTER FITTING" << std::endl;
			std::cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<std::endl;
			std::cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<std::endl;
			std::cout << "[INFO] NLL + c = " <<  minNll << std::endl;
			std::cout << "-----------------------" << std::endl;
		}
			
		if (minNll < global_minNll){
      global_minNll = minNll;
			snap.assignValueOnly(*params);
      best_index=id;
		}
	}
  cat->setIndex(best_index);
	params->assignValueOnly(snap);
	
	if (!silent) {
		std::cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<std::endl;
		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}

FTestCategoryResult runEnvelopeFTest(const FTestCategoryInputs& input, BackgroundFitOptions& options) {
  FTestCategoryResult result;
		// Standard F-Test to find the truth functions
			for (vector<string>::const_iterator funcType=input.functionClasses.begin(); funcType!=input.functionClasses.end(); funcType++){

			double thisNll=0.; double prevNll=0.; double chi2=0.; double prob=0.; 
			int order=1; int prev_order=0; int cache_order=0;

			RooAbsPdf *prev_pdf=NULL;
			RooAbsPdf *cache_pdf=NULL;
			std::vector<int> pdforders;

      string typePrefix;
      // while (prob<options.upperEnvThreshold){
			// while (prob<options.upperEnvThreshold && order < 7){
      //   RooAbsPdf *bkgPdf = getPdf(*input.pdfsModel,*funcType,order,&typePrefix,Form("ftest_pdf_%d_%s",(cat+input.catOffset),input.energyExt.c_str()), input.catName, /*sideband*/false);
      //   cout << "[INFO] Try to find " << *funcType << " " << order << " " << bkgPdf << endl;
      //   // if (order < 5){ //DEBUG
      //   //   order++;
      //   //   continue;
      //   // } 
      //   // bkgPdf->Print("v");bkgPdf_SB->Print("v");
      //   if (!bkgPdf){
      //     // assume this order is not allowed
      //     cout << "[INFO] bkg pdf is not exist " << endl;
      //     order++;
      //   }
      //   else {
      //     cout << "[INFO] get pdf called " << *funcType << " " << order << " " << bkgPdf << endl;
      //     int fitStatus = 0;
      //     bkgPdf->Print();
      //     bkgPdf_SB->Print();
      //     // runFit(bkgPdf,bkgPdf_SB,input.data,input.dataSB,&thisNll,&fitStatus,/*max iterations*/3, /*options.blindFit*/options.blindFit, options.flashggCats[cat].c_str(), typePrefix, order, options.nBinsForMass, bkgPdf);
      //     // if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<std::endl;       
      //     // chi2 = 2.*(prevNll-thisNll);
      //     // if (chi2<0. && order>1) chi2=0.;
      //     // if (prev_pdf!=NULL){
      //     //   prob = getProbabilityFtest(chi2,order-prev_order,prev_pdf,bkgPdf,mass,input.data
      //     //       ,Form("%s/Ftest_from_%s%d_cat%d.pdf",options.outDir.c_str(),funcType->c_str(),order,(cat+input.catOffset)));
      //     //   std::cout << "[INFO]  F-test Prob(chi2>chi2(input.data)) == " << prob << std::endl;
      //     // } else {
      //     //   prob = 0;
      //     // }
      //     // double gofProb=0;
      //     // // otherwise we get it later ...
      //     // if (!options.saveMultiPdf) plot(mass,bkgPdf,input.data,Form("%s/%s%d_cat%d.png",options.outDir.c_str(),funcType->c_str(),order,(cat+input.catOffset)),options.flashggCats,fitStatus,&gofProb);
      //     // cout << "[INFO]\t funcType prev_order order prevNll thisNll chi2 prob" << endl;
      //     // cout << "[INFO]\t " << *funcType << " " << prev_order << " " << order << " " << prevNll << " " << thisNll << " " << chi2 << " " << prob << endl;
      //     // //fprintf(input.resFile,"%15s && %d && %10.2f && %10.2f && %10.2f \\\\\n",funcType->c_str(),order,thisNll,chi2,prob);
      //     // prevNll=thisNll;
      //     // if (prev_pdf==NULL) {
      //     //   cache_pdf=bkgPdf;
      //     //   cache_order=order;
      //     // }
      //     // else {
      //     //   cache_pdf=prev_pdf;
      //     //   cache_order=prev_order;
      //     // }
      //     // prev_order=order;
      //     // prev_pdf=bkgPdf;
      //     // cout << "[INFO] Ftest\t " << *funcType << " " << cache_order << " " << prev_order << endl;
      //     order++;
      //   }
      // }

			// fprintf(input.resFile,"%15s & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
      // cout << "[INFO] Ftest upper limit " << cache_order << " " << prob << endl;
      // if (cache_order > 0) {
			//   result.choices.insert(pair<string,int>(*funcType,cache_order));
      //   result.pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
      // }
			// int truthOrder = cache_order;

			// Now run loop to determine functions inside envelope
			if (options.saveMultiPdf || options.objectiveOnly){
				chi2=0.;
				thisNll=0.;
				prevNll=0.;
				prob=0.;
				order=options.minPdfOrder;
				prev_order=0;
				cache_order=0;
				prev_pdf=NULL;
				cache_pdf=NULL;
					std::cout << "[INFO] Determining Envelope Functions for Family " << *funcType << ", cat " << input.cat << std::endl;
				std::cout << "[INFO] Upper end Threshold for highest order function " << options.upperEnvThreshold <<std::endl;

        while (order <= options.maxPdfOrder){ cout << "[DEBUG INFO] Try all functions " << endl;// DEBUG
						if (!options.candidateOrders.empty() &&
						    std::find(options.candidateOrders.begin(), options.candidateOrders.end(), order) == options.candidateOrders.end()) {
						  order++;
						  continue;
						}
						RooAbsPdf *bkgPdf = getPdf(*input.pdfsModel,*funcType,order,&typePrefix,Form("env_pdf_cat%d_%s",(input.cat+input.catOffset),input.energyExt.c_str()), input.catName, /*sideband*/false, options.useDiGaussian);
          cout << "[INFO] Try to find " << *funcType << " " << order << " " << bkgPdf << endl;
          // if (order != 3){ cout << "[DEBUG INFO] Just look into order 2 " << endl;// DEBUG
          //   order++;
          //   continue;
          // } 
          
          if (!bkgPdf){
            // assume this order is not allowed
            cout << "[INFO] bkg pdf is not exist " << endl;
            order++;
            continue;
          }
          // bkgPdf->Print("v");bkgPdf_SB->Print("v");
          fprintf(input.dfile,"\n%s %d\n",funcType->c_str(),order);
          cout << "[INFO] get pdf called " << *funcType << " " << order << " " << bkgPdf << endl;
          //RooFitResult *fitRes;
          int fitStatus=0;
	          int covQual = -999;
	          double edm = std::numeric_limits<double>::infinity();
	          int boundaryHit = -1;
	          runFit(bkgPdf, input.data, input.dataSB, &thisNll,
	                 &fitStatus, options, options.blindFit,
	                 options.flashggCats[input.cat].c_str(), typePrefix, order,
	                 options.nBinsForMass, bkgPdf, options.doChi2Fit,
	                 &covQual, &edm, &boundaryHit);
          //thisNll = fitRes->minNll();
          const std::string fitLabel = std::string(bkgPdf->GetName()) + "_SB";
          if (fitStatus!=0) std::cout << "[WARNING] Warning -- Fit status for " << fitLabel.c_str() << " at " << fitStatus <<std::endl;
	          std::cout << "[CANDIDATE_FIT] category=" << input.catName
	                    << " family=" << *funcType
	                    << " order=" << order
	                    << " pdf=" << fitLabel.c_str()
	                    << " status=" << fitStatus
	                    << " cov_qual=" << covQual
	                    << " edm=" << edm
	                    << " min_nll=" << thisNll
	                    << " boundary_hit=" << boundaryHit << std::endl;
	          // Keep the fit-quality result available independently of whether
	          // the active campaign uses it as a candidate-selection gate.
	          const bool diagnosticFitQualityPass =
	                                        options.passesNumericalFitQuality(fitStatus, covQual, edm);
	          const bool passedFitQuality =
	                                        !options.gateFitQuality || diagnosticFitQualityPass;
	          std::cout << "[CANDIDATE_DIAGNOSTIC] category=" << input.catName
	                    << " family=" << *funcType << " order=" << order
	                    << " stage=fit_quality"
	                    << " fit_status=" << fitStatus
	                    << " cov_qual=" << covQual
	                    << " edm=" << edm
	                    << " boundary_hit=" << boundaryHit
	                    << " boundary_gating=0"
	                    << " fit_status_gating=" << (options.requireFitStatusZero ? 1 : 0)
	                    << " fit_quality_gating=" << (options.gateFitQuality ? 1 : 0)
	                    << " diagnostic_passed=" << (diagnosticFitQualityPass ? 1 : 0)
	                    << " passed=" << (passedFitQuality ? 1 : 0) << std::endl;
	          if (!passedFitQuality) {
	            std::cout << "[CANDIDATE_REJECT] category=" << input.catName
	                      << " family=" << *funcType << " order=" << order
	                      << " stage=fit_quality accepted=0" << std::endl;
	            order++;
	            continue;
	          }
          double myNll = 2.*thisNll;
          chi2 = max(2.*(prevNll-thisNll), 0.);
          // if (chi2<0. && order>1) chi2=0.;
          if (prev_order == 0) prob = 0; // lowest order do not need F-test
          else prob = TMath::Prob(chi2,order-prev_order);

          cout << "[INFO] \t " << fitLabel.c_str() << ", prevNll = " << prevNll << ", thisNll = " << thisNll << ", chi2 = " << chi2 << ", prob = " << prob << endl;
          fprintf(input.dfile,"%s , order = %d, prev_order = %d, prevNll = %f, thisNll = %f, prob = %f\n",fitLabel.c_str(),order,prev_order,prevNll,thisNll,prob);
          if (prev_pdf==NULL) {
            cache_pdf=bkgPdf;
            cache_order=order;
          }
          else {
            cache_pdf=prev_pdf;
            cache_order=prev_order;
          }
          
          // Calculate goodness of fit for the thing to be included (will use toys for lowstats)!
          double gofProb =0; 
          // if(!options.blindFit) plot(mass,/*RooAbsPdf*/bkgPdf,input.data,Form("%s/%s%d_cat%d.png",options.outDir.c_str(),funcType->c_str(),order,(cat+input.catOffset)),options.flashggCats,fitStatus,&gofProb);
          // else plot(mass,/*RooAbsPdf*/bkgPdf_SB,input.dataSB,Form("%s/%s%d_cat%d.png",options.outDir.c_str(),funcType->c_str(),order,(cat+input.catOffset)),options.flashggCats,fitStatus,&gofProb);
	          plot(input.mass, /*RooAbsPdf*/bkgPdf, input.data, Form("%s/%s_%s%d.png", options.outDir.c_str(), input.catName.c_str(), funcType->c_str(), order), options.flashggCats, fitStatus, &gofProb, options);
          // RooChi2Var chi2Var(Form("chi2Var_%s",bkgPdf->GetName()), Form("chi2Var_%s",bkgPdf->GetName()), *bkgPdf_SB, *thisdataSBBinned);
          // double chi2 = chi2Var.getVal();

          fprintf(input.dfile,"p-value(chi2) = %.2f, ",gofProb);
	          std::cout << "[CANDIDATE_GOF] category=" << input.catName
	                    << " family=" << *funcType << " order=" << order
	                    << " gof_prob=" << gofProb << std::endl;
	          if (gofProb <= options.minimumGofProbability) {
            std::cout << "[INFO] " << *funcType << " with order " << order
                      << " failed GOF threshold " << options.minimumGofProbability
                      << ", skipping KS/SST/F-test" << std::endl;
            order++;
            continue;
          }
          // plot(mass,/*RooAbsPdf*/bkgPdf_SB,input.dataSB,Form("%s/%s%dSB_cat%d.png",options.outDir.c_str(),funcType->c_str(),order,(cat+input.catOffset)),options.flashggCats,fitStatus,&gofProb);
                      
          // Calculate KS test probability
	          double ksProb = getKSProb(input.mass, bkgPdf, input.dataFull, Form("%s/%s_%s%d", options.outDir.c_str(), input.catName.c_str(), funcType->c_str(), order), options);
          cout << "[INFO] \t KS test probability = " << ksProb << endl;
	          std::cout << "[CANDIDATE_KS] category=" << input.catName
	                    << " family=" << *funcType << " order=" << order
	                    << " ks_prob=" << ksProb << std::endl;
          fprintf(input.dfile,"KS test probability = %.2f, ",ksProb);
	          if (ksProb <= options.minimumKsProbability) {
            std::cout << "[INFO] " << *funcType << " with order " << order
                      << " failed KS threshold " << options.minimumKsProbability
                      << ", skipping SST/F-test" << std::endl;
            order++;
            continue;
          }

          bool passedSSTest = true;
          if (options.skipSSTest) {
            std::cout << "[INFO] Skipping SpuriousSignalTest for "
                      << input.catName << " " << *funcType << order << std::endl;
          } else {
            try {
              passedSSTest = SpuriousSignalTest(input.mass, TString(input.catName.c_str()),
                                                TString(funcType->c_str()), order,
                                                options.nsigIn, input.runPeriod, options);
            } catch (const std::exception& e) {
              std::cerr << "[ERROR] SpuriousSignalTest failed for "
                        << input.catName << " " << *funcType << order
                        << ": " << e.what() << std::endl;
              passedSSTest = false;
            }
          }
          fprintf(input.dfile,"sstest status = %d, ",passedSSTest);
	          std::cout << "[CANDIDATE_SST] category=" << input.catName
	                    << " family=" << *funcType << " order=" << order
	                    << " skipped=" << (options.skipSSTest ? 1 : 0)
	                    << " passed=" << (passedSSTest ? 1 : 0) << std::endl;
          if (!passedSSTest) {
            std::cout << "[INFO] " << *funcType << " with order " << order
                      << " failed SpuriousSignalTest, skipping F-test" << std::endl;
            order++;
            continue;
          }

          const bool passedFTest = prob < options.upperEnvThreshold;
          std::cout << "[CANDIDATE_FTEST] category=" << input.catName
                    << " family=" << *funcType << " order=" << order
                    << " ftest_prob=" << prob
                    << " threshold=" << options.upperEnvThreshold
                    << " passed=" << (passedFTest ? 1 : 0) << std::endl;
          if (passedFTest) { // Looser requirements for the envelope
          // if (1) { cout << "[DEBUG INFO] Ignore F-test " << endl;// DEBUG

            // if (gofProb > 0.01 || order == truthOrder ) {  // Good looking fit or one of our regular truth functions
            if (1) {  // GOF, KS, and SST have already passed in sequence.
            // if (1) { cout << "[DEBUG INFO] Ignore chi2 and KS test "// DEBUG

              std::cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " "<< gofProb 
                << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() << std::endl;
              result.allPdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),order),bkgPdf));
              result.storedPdfs.add(*bkgPdf);
              pdforders.push_back(order);
              prev_order=order;
              prev_pdf=bkgPdf;
              prevNll=thisNll;
              cout << "[INFO] Ftest envelope " << *funcType << " " << order << " " << prob << endl;
	              std::cout << "[CANDIDATE_SELECT] category=" << input.catName
	                        << " family=" << *funcType << " order=" << order
	                        << " ftest_prob=" << prob << " accepted=1" << std::endl;
              fprintf(input.dfile,"Ftest envelope = %.2f",prob);
          
              // Keep track but we shall redo this later
              if ((myNll + bkgPdf->getVariables()->getSize()) < result.minimumNLLSoFar) {
                result.simpleBestFitPdfIndex = result.storedPdfs.getSize()-1;
                result.minimumNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
              }
            }

					}
          order++;
				}

        if (cache_order > 0) {
          result.choices.insert(pair<string,int>(*funcType,cache_order));
          result.pdfs.insert(pair<string,RooAbsPdf*>(Form("%s%d",funcType->c_str(),cache_order),cache_pdf));
          if (input.resFile && cache_order == options.maxPdfOrder){
            fprintf(input.resFile,"%15s F-test & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order,chi2,prob);
          } else if (input.resFile) {
            fprintf(input.resFile,"%15s Fail F-test & %d & %5.2f & %5.2f \\\\\n",funcType->c_str(),cache_order+1,chi2,prob);
          }
        }
        result.choicesEnvelope.insert(pair<string,std::vector<int> >(*funcType,pdforders));
			}
		}

			if (input.resFile) fprintf(input.resFile,"\\hline\n");
	    for (map<string,RooAbsPdf*>::iterator it=result.pdfs.begin(); it!=result.pdfs.end(); it++){
	      std::cout << "[INFO]  " << it->first << std::endl;
	    }
	    // Some categories don't have any functions in the envelope
			plot(input.mass, /*map<string,RooAbsPdf*>*/result.pdfs, input.data, Form("%s/truths_%s", options.outDir.c_str(), input.catName.c_str()), options.flashggCats, input.cat, -1, options);

  return result;
}
