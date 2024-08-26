#ifndef WANALYZER_H
#define WANALYZER_H

#include "wPair.h"
#include "wFile.h"
#include "wMode.h"
#include "wBMP.h"
#include "wContour.h"

namespace wAnalyzer{
    void curvatureAnalysisFromImageJExtractedOrgan(std::string,std::string,int,int,int,bool);
}
#endif