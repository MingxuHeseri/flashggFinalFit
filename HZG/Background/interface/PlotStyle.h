#pragma once

#include "TString.h"
class TPad;

void setTDRStyle();
void CMS_lumi(TPad*, int, int, TString);
extern bool writeExtraText;
extern TString extraText, lumi_8TeV, lumi_7TeV, lumi_sqrtS;
