#include "../interface/PdfFactory.h"

#include <cstdlib>
#include <iostream>
#include <vector>

#include "boost/property_tree/json_parser.hpp"
#include "boost/property_tree/ptree.hpp"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "TIterator.h"
#include "TString.h"

using namespace std;

// Save parameters to JSON file
bool saveParamsToJSON(const string& cat, const string& type, int order, const RooAbsPdf* pdf) {
  // Create params directory (if it doesn't exist)
  system("mkdir -p params");
  
  std::string fileName = Form("params_%s_%s%d.json", cat.c_str(), type.c_str(), order);
  std::string filePath = "params/" + fileName;
  
  try {
    // Get PDF parameters
    RooArgSet* params = pdf->getParameters(RooArgSet());
    boost::property_tree::ptree pt;
    
    // Iterate through parameters and add to JSON tree
    for (RooAbsArg *arg : *params) {
      RooRealVar *param = dynamic_cast<RooRealVar*>(arg);
      if (!param) continue;
      if (!param->isConstant()) {  // Only save floating parameters
        pt.put(boost::property_tree::ptree::path_type(param->GetName(), '\0'), param->getVal());
      }
    }
    delete params;
    
    // Write to JSON file
    boost::property_tree::write_json(filePath, pt);
    std::cout << "[INFO] Saved parameters to " << filePath << std::endl;
    return true;
  }
  catch (const std::exception& e) {
    std::cerr << "[ERROR] Failed to save parameters to " << filePath << ": " << e.what() << std::endl;
    return false;
  }
}

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, string *typePrefix, const char* ext, const string& cat, bool doBlind, bool useDiGaussian){
  RooAbsPdf* pdf = nullptr;
  pdfsModel.setParameterConfigCategory(cat);
  
  // First create the PDF instance
  if (type=="Bernstein") {
    pdf = pdfsModel.getBernstein(Form("%s_bern%d",ext,order),order,doBlind,
                                 "Bernstein");
    *typePrefix = "bern"; 
  }
  else if (type=="BernsteinPositive") {
    pdf = pdfsModel.getBernstein(Form("%s_bernpositive%d",ext,order),order,doBlind,
                                 "BernsteinPositive");
    *typePrefix = "bernpositive";
  }
  else if (type=="Chebychev") {
    pdf = pdfsModel.getChebychev(Form("%s_cheb%d",ext,order),order);
    *typePrefix = "cheb";
  }
  else if (type=="Exponential") {
    pdf = pdfsModel.getExponential(Form("%s_exp%d",ext,order),order, doBlind);
    *typePrefix = "exp";
  }
  else if (type=="PowerLaw") {
    pdf = pdfsModel.getPowerLaw(Form("%s_pow%d",ext,order),order, doBlind);
    *typePrefix = "pow";
  }
  else if (type=="Laurent") {
    pdf = pdfsModel.getLaurentSeries(Form("%s_lau%d",ext,order),order);
    *typePrefix = "lau";
  }
  else if (type=="BernsteinStepxGauFast") {
    if (order < 2) return NULL;
    pdf = pdfsModel.getBernsteinStepxGauFast(Form("%s_bernfast%d",ext,order-1),order, doBlind, useDiGaussian);
    *typePrefix = useDiGaussian ? "bernfastdigau" : "bernfast";
  }
  else if (type=="ExpPolyTail" || type=="PowerLogTail" || type=="LaurentTail") {
    if (order < 1 || order > 3 || useDiGaussian) return NULL;
    *typePrefix = type=="ExpPolyTail" ? "exppolytail" :
                  (type=="PowerLogTail" ? "powerlogtail" : "laurenttail");
    pdf = pdfsModel.getEmpiricalTail(
        Form("%s_%s%d", ext, typePrefix->c_str(), order), order, doBlind, type);
  }
  else if (type=="BernsteinErfTurnOn") {
    if (order < 2 || order > 6 || useDiGaussian) return NULL;
    pdf = pdfsModel.getBernsteinErfTurnOn(
        Form("%s_bernerf%d", ext, order-1), order, doBlind);
    *typePrefix = "bernerf";
  }
  else if (type=="ExponentialErfTurnOn") {
    if (order < 1 || order > 5 || order % 2 == 0 || useDiGaussian) return NULL;
    pdf = pdfsModel.getExponentialErfTurnOn(
        Form("%s_experf%d", ext, order/2+1), order, doBlind);
    *typePrefix = "experf";
  }
  else if (type=="PowerLawErfTurnOn") {
    if (order < 1 || order > 5 || order % 2 == 0 || useDiGaussian) return NULL;
    pdf = pdfsModel.getPowerLawErfTurnOn(
        Form("%s_powerf%d", ext, order/2+1), order, doBlind);
    *typePrefix = "powerf";
  }
  else if (type=="LaurentErfTurnOn") {
    if (order < 2 || order > 6 || useDiGaussian) return NULL;
    pdf = pdfsModel.getLaurentErfTurnOn(
        Form("%s_lauerf%d", ext, order), order, doBlind);
    *typePrefix = "lauerf";
  }
  else if (type=="ExpBernsteinErfTurnOn") {
    if (order < 1 || order > 6 || useDiGaussian) return NULL;
    pdf = pdfsModel.getExpBernsteinErfTurnOn(
        Form("%s_expbernerf%d", ext, order), order, doBlind);
    *typePrefix = "expbernerf";
  }
  else if (type=="ExpChebychevErfTurnOn") {
    if (order < 1 || order > 8 || useDiGaussian) return NULL;
    pdf = pdfsModel.getExpChebychevErfTurnOn(
        Form("%s_expcheberf%d", ext, order), order, doBlind);
    *typePrefix = "expcheberf";
  }
  else if (type=="BernsteinStepxGauFastUnconstrained") {
    if (order < 2 || order > 6 || useDiGaussian) return NULL;
    pdf = pdfsModel.getBernsteinStepxGauFastUnconstrained(
        Form("%s_bernfastunconstrained%d", ext, order - 1), order, doBlind);
    *typePrefix = "bernfastunconstrained";
  }
  else if (type=="BernsteinStepxGauAnalytic") {
    if (order < 2 || order > 6 || useDiGaussian) return NULL;
    pdf = pdfsModel.getBernsteinStepxGauAnalytic(
        Form("%s_bernanalytic%d", ext, order-1), order, doBlind);
    *typePrefix = "bernanalytic";
  }
  else if (type=="TurnOnExpChebychev") {
    if (order < 1 || order > 8 || useDiGaussian) return NULL;
    pdf = pdfsModel.getTurnOnExpChebychev(
        Form("%s_turnonexpcheb%d", ext, order), order, doBlind);
    *typePrefix = "turnonexpcheb";
  }
  else if (type=="HZGErfExp") {
    if (order != 1 || useDiGaussian) return NULL;
    pdf = pdfsModel.getHZGErfExp(Form("%s_hzgerfexp", ext), doBlind);
    *typePrefix = "hzgerfexp";
  }
  else if (type=="HZGPowerExp") {
    if (order != 1 || useDiGaussian) return NULL;
    pdf = pdfsModel.getHZGPowerExp(Form("%s_hzgpowerexp", ext), doBlind);
    *typePrefix = "hzgpowerexp";
  }
  else if (type=="HZGErfBernstein") {
    if (order < 2 || order > 3 || useDiGaussian) return NULL;
    pdf = pdfsModel.getHZGErfBernstein(
        Form("%s_hzgerfbern%d", ext, order), order, doBlind);
    *typePrefix = "hzgerfbern";
  }
  else if (type=="Fkp") {
    if (order < 0 || order > 2 || useDiGaussian) return NULL;
    pdf = pdfsModel.getFkp(Form("%s_fkp%d", ext, order), order, doBlind);
    *typePrefix = "fkp";
  }
  else if (type=="ExponentialStepxGau") {
    pdf = pdfsModel.getExponentialStepxGau(Form("%s_exp%d",ext,order/2+1),order, doBlind, useDiGaussian, /*const_cp1*/true);
    *typePrefix = useDiGaussian ? "expdigau" : "exp";
  }
  else if (type=="PowerLawStepxGau") {
    pdf = pdfsModel.getPowerLawStepxGau(Form("%s_pow%d",ext,order/2+1),order, doBlind, useDiGaussian, /*const_cp1*/true);
    *typePrefix = useDiGaussian ? "powdigau" : "pow";
  }
  else if (type=="LaurentStepxGau") {
    // Keep the candidate family aligned with the Run 3 background scan:
    // LaurentStepxGau starts at order 2 and is scanned through order 6.
    if (order < 2 || order > 6) return NULL;
    // 根据cat的内容设置p值数组
    std::vector<int> ps;
    if (cat.find("ggF1") != std::string::npos) {
        ps = {-2, -9, -11, -12};
    } else if (cat.find("ggF2") != std::string::npos) {
        ps = {-5, -11, -12, -13};
    } else if (cat.find("ggF3") != std::string::npos) {
        ps = {-7, -8, -5, -9};
    } else if (cat.find("ggF4") != std::string::npos) {
        ps = {-4, -8, -5, -9};
    } else if (cat.find("VBF0") != std::string::npos) {
        ps = {-3, -4, -6, -5};
    } else if (cat.find("VBF1") != std::string::npos) {
        ps = {-11, -4, -5, -12};
    } else if (cat.find("VBF2") != std::string::npos) {
        ps = {-4, -13, -15, -5};
    } else if (cat.find("VBF3") != std::string::npos) {
        ps = {-2, -9, -10, -11};
    } else {  // 默认情况，保持原来的逻辑
        ps = {-4, -5, -6, -7};
    }
    pdf = pdfsModel.getLaurentStepxGau(Form("%s_lau%d",ext,order),order, doBlind, useDiGaussian, /*const_cp1*/true);
    *typePrefix = useDiGaussian ? "laudigau" : "lau";
  }
  else if (type=="ModGaus"){
    if (order == 1) {
      pdf = pdfsModel.getModGaus(Form("%s_modgaus",ext), doBlind);
      *typePrefix = "modgaus";
    }
    else return NULL;
  }
  else if (type=="ExpModGauss"){
    if (order == 1) {
      pdf = pdfsModel.getExpModGaussian(Form("%s_modgauss",ext), doBlind);
      *typePrefix = "modgauss";
    }
    else return NULL;
  }
  else if (type=="AsymGenGauss"){
    if (order == 1) {
      pdf = pdfsModel.getAsymGenGaussian(Form("%s_asymgauss",ext), doBlind);
      *typePrefix = "asymgauss";
    }
    else return NULL;
  }
  else {
    cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
    return NULL;
  }
  
  return pdf;
}

std::map<std::string, std::string> defaultNamingMap() {
  std::map<std::string, std::string> namingMap;
  namingMap.insert(std::pair<std::string, std::string>("Bernstein", "pol"));
  namingMap.insert(std::pair<std::string, std::string>("BernsteinPositive", "bernpositive"));
  namingMap.insert(std::pair<std::string, std::string>("Exponential", "exp"));
  namingMap.insert(std::pair<std::string, std::string>("ExpPolyTail", "exppolytail"));
  namingMap.insert(std::pair<std::string, std::string>("PowerLogTail", "powerlogtail"));
  namingMap.insert(std::pair<std::string, std::string>("LaurentTail", "laurenttail"));

  namingMap.insert(std::pair<std::string, std::string>("PowerLaw", "pow"));
  namingMap.insert(std::pair<std::string, std::string>("Laurent", "lau"));
  namingMap.insert(std::pair<std::string, std::string>("BernsteinStepxGauFast", "bernfast"));
  namingMap.insert(std::pair<std::string, std::string>("BernsteinErfTurnOn", "bernerf"));
  namingMap.insert(std::pair<std::string, std::string>("ExponentialErfTurnOn", "experf"));
  namingMap.insert(std::pair<std::string, std::string>("PowerLawErfTurnOn", "powerf"));
  namingMap.insert(std::pair<std::string, std::string>("LaurentErfTurnOn", "lauerf"));
  namingMap.insert(std::pair<std::string, std::string>("ExpBernsteinErfTurnOn", "expbernerf"));
  namingMap.insert(std::pair<std::string, std::string>("ExpChebychevErfTurnOn", "expcheberf"));
  namingMap.insert(std::pair<std::string, std::string>("BernsteinStepxGauFastUnconstrained", "bernfastunconstrained"));
  namingMap.insert(std::pair<std::string, std::string>("BernsteinStepxGauAnalytic", "bernanalytic"));
  namingMap.insert(std::pair<std::string, std::string>("TurnOnExpChebychev", "turnonexpcheb"));
  namingMap.insert(std::pair<std::string, std::string>("HZGErfExp", "hzgerfexp"));
  namingMap.insert(std::pair<std::string, std::string>("HZGPowerExp", "hzgpowerexp"));
  namingMap.insert(std::pair<std::string, std::string>("HZGErfBernstein", "hzgerfbern"));
  namingMap.insert(std::pair<std::string, std::string>("Fkp", "fkp"));
  namingMap.insert(std::pair<std::string, std::string>("ExponentialStepxGau", "exp"));
  namingMap.insert(std::pair<std::string, std::string>("PowerLawStepxGau", "pow"));
  namingMap.insert(std::pair<std::string, std::string>("LaurentStepxGau", "lau"));
  namingMap.insert(std::pair<std::string, std::string>("ExpModGauss", "modgau"));
  namingMap.insert(std::pair<std::string, std::string>("AsymGenGauss", "asymgau"));
  return namingMap;
}
