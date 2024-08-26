/*************************************************************************/
// 2D Organ contour analysis
//Author: Zining Wang (wangzining16@mails.ucas.ac.cn/wangzining2020@g.ecc.u-tokyo.ac.jp)
//To read, analyze, and plot the txt of organ contour extracted from ImageJ 
/*********************************************************************** */
#include "../include/wFile.h"
#include "../include/wContour.h"
#include "../include/wMode.h"
#include "../include/wBMP.h"
#include "../include/wAnalyzer.h"

std::string inputFileDirectory="/mnt/c/again_and_again/processed_experiment_data/2024_07_24_oxalis_debilis_triangular/contour/";
std::string outputFileDirectory="/mnt/c/again_and_again/processed_experiment_data/2024_07_24_oxalis_debilis_triangular/postprocess/";
std::string modeFile = "../mode.txt";

int main(){
    wMode wmode;
    wmode.readFromFile(modeFile);
    
    if(wmode.get_major_mode()=="trial"){
    }
    else if(wmode.get_major_mode()=="curvature_analysis"){
        if(wmode.get_minor_mode()=="single"){
            int boundaryPointNum = 100;
            int pointsAway = 10;
            int numAveraging = 3;
            bool debugMode = 0;
            std::string fileName = "OD-1-1";
            std::string fileFormat = ".csv";
            wAnalyzer::curvatureAnalysisFromImageJExtractedOrgan(inputFileDirectory+fileName+fileFormat,fileName,boundaryPointNum,pointsAway,numAveraging,debugMode);
        }
        
        else if(wmode.get_minor_mode()=="multiple"){
            int boundaryPointNum = 100;
            int pointsAway = 10;
            int numAveraging = 3;
            bool debugMode = 0;
            std::vector<std::string> fileNames = {"OD-1-1","OD-1-3"};
            std::string fileFormat = ".csv";
            for(int i=0;i<fileNames.size();i++){
                std::string fileName = fileNames[i];
                wAnalyzer::curvatureAnalysisFromImageJExtractedOrgan(inputFileDirectory+fileName+fileFormat,fileName,boundaryPointNum,pointsAway,numAveraging,debugMode);
            }
        }
        else if(wmode.get_minor_mode()=="batch"){
            int boundaryPointNum = 100;
            int pointsAway = 10;
            int numAveraging = 3;
            bool debugMode = 0;
            std::vector<std::string> fileNames = wFile::getFilesNameInDirectory(inputFileDirectory);
            std::string fileFormat = ".csv";
            for(int i=0;i<fileNames.size();i++){
                std::string fileName = fileNames[i];
                std::cout<<"Read and analyze curvature from "<<inputFileDirectory+fileName+fileFormat<<std::endl;
                wAnalyzer::curvatureAnalysisFromImageJExtractedOrgan(inputFileDirectory+fileName+fileFormat,fileName,boundaryPointNum,pointsAway,numAveraging,debugMode);
            }
        }
        else{
            std::cerr<<"No mode selected for minor_mode ("<<wmode.get_minor_mode()<<")"<<std::endl;
            exit(1);
        }
    }
    else{
        std::cerr<<"No mode selected for major_mode ("<<wmode.get_major_mode()<<")"<<std::endl;
        exit(1);
    }


    wmode.recordToFile(outputFileDirectory+"mode.txt");
    return 0;
}