#ifndef FLASHGGFINALFIT_HZG_PDF_FACTORY_H
#define FLASHGGFINALFIT_HZG_PDF_FACTORY_H

#include <map>
#include <string>

#include "RooAbsPdf.h"
#include "PdfModelBuilder.h"

bool saveParamsToJSON(const std::string& cat, const std::string& type, int order, const RooAbsPdf* pdf);
RooAbsPdf* getPdf(PdfModelBuilder& pdfsModel, std::string type, int order, std::string* typePrefix, const char* ext="", const std::string& cat="", bool doBlind=false, bool useDiGaussian=false);
std::map<std::string, std::string> defaultNamingMap();

#endif
