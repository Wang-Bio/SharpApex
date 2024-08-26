#include "../include/wAnalyzer.h"

namespace wAnalyzer{
    void curvatureAnalysisFromImageJExtractedOrgan(std::string inputFile, std::string fileName,int boundaryPointNum,int pointsAway,int numAveraging,bool debugMode){
        //1. read and draw raw contour
        std::vector<std::pair<double,double>> rawContour = wFile::readContour_raw(inputFile);
            if(debugMode==1){wFile::outputContour(rawContour,outputFileDirectory+fileName+"_raw_contour.txt");};
            //if(debugMode==1){wBMP::drawContour(rawContour,outputFileDirectory+fileName+"_raw_contour.bmp");};
        //2. raw contour in imageJ need to be flipped
        std::vector<std::pair<double,double>> flipContour = wContour::flipY(rawContour);
            wBMP::drawContour(flipContour,outputFileDirectory+fileName+"_raw_contour.bmp");
        //3. raw contour need to be normalized
        std::vector<std::pair<double,double>> normContour = wContour::normalization(flipContour);
            if(debugMode==1){wFile::outputContour(normContour,outputFileDirectory+fileName+"norm_contour.txt");};
            //if(debugMode==1){wBMP::drawContour(normContour,outputFileDirectory+fileName+"_norm_contour.bmp");};
        //4. normalized contour could be parameterized by eg., 100 points with equal distance
        std::vector<std::pair<double,double>> paramContour = wContour::parameterization(normContour,boundaryPointNum);
            wFile::outputContour(paramContour,outputFileDirectory+fileName+"param_contour.txt");
            //if(debugMode==1){wBMP::drawContour(paramContour,outputFileDirectory+fileName+"_param_contour.bmp");};
        //5. curvature could be calculated from the parameterized contour
        std::vector<double> paramCurvature = wContour::vec_curvature_circleFittingKasa_threeBoundaryPoints(paramContour,pointsAway,numAveraging);
            wFile::vd_fout(paramCurvature,outputFileDirectory+fileName+"_curvature.txt");
            wBMP::drawCurvature(paramCurvature,outputFileDirectory+fileName+"_curvature.bmp");
    }
}