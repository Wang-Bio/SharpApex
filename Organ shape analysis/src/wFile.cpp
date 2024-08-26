#include "../include/wFile.h"

namespace wFile{
    std::vector<std::pair<double,double>> readVecPairDouble(const std::string& inputFile){
                bool debugMode=0;
                if(debugMode==1){std::cout<<"debugMode for readVecPairDouble"<<std::endl;};
        std::vector<std::pair<double,double>> PairDouble;
        std::ifstream fin(inputFile,std::ios::in);
        if(!fin.is_open()){
            std::cerr<<"Error: the inputFile not found ("<<inputFile<<")"<<std::endl;
            exit(1);
        }
        while(!fin.eof()){
            double x_tmp,y_tmp;
            fin>>x_tmp;
            fin>>y_tmp;
            PairDouble.push_back(std::make_pair(x_tmp,y_tmp));
                if(debugMode==1){std::cout<<PairDouble.size()<<" "<<x_tmp<<","<<y_tmp<<std::endl;};
        }

        return PairDouble;
    }


    std::vector<std::pair<double,double>> readContour_posi_only(const std::string& inputFile){
        std::vector<std::pair<double,double>> rawContour;
            bool debugMode=0;
        std::ifstream fin(inputFile,std::ios::in);
        if(!fin.is_open()){
            std::cerr<<"Error: the inputFile not found ("<<inputFile<<")"<<std::endl;
            exit(1);
        }
        std::string headerLine;
        std::getline(fin,headerLine); //headerLine is should be discarded
                if(debugMode==1){std::cout<<"HeaderLine "<<headerLine<<std::endl;};
        while(!fin.eof()){
            int posi_x,posi_y;
            fin>>posi_x;
            fin>>posi_y;
            std::pair<double,double> contourPointPosition = std::make_pair(posi_x,posi_y);
            rawContour.push_back(contourPointPosition);
                if(debugMode==1){std::cout<<"for point "<<rawContour.size()<<" its posi is "<<contourPointPosition.first<<","<<contourPointPosition.second<<std::endl;};
        }
                if(debugMode==1){std::cout<<"File "<<inputFile<<" has "<<rawContour.size()<<" contour points"<<std::endl;};
        return rawContour;
    }

    std::vector<std::pair<double,double>> readContour_raw(const std::string& inputFile){
        std::vector<std::pair<double,double>> rawContour;
            bool debugMode=0;
        std::ifstream fin(inputFile,std::ios::in);
        if(!fin.is_open()){
            std::cerr<<"Error: the inputFile not found ("<<inputFile<<")"<<std::endl;
            exit(1);
        }
        std::string headerLine;
        std::getline(fin,headerLine); //headerLine is should be discarded
                if(debugMode==1){std::cout<<"HeaderLine "<<headerLine<<std::endl;};
        std::string line;
        while(std::getline(fin,line)){
            std::stringstream ss(line);
            std::string token;
            int index;
            double posi_x, posi_y;
            std::getline(ss,token,',');
            index = std::stod(token);
            std::getline(ss,token,',');
            posi_x = std::stod(token);
            std::getline(ss,token,',');
            posi_y = std::stod(token);
            std::pair<double,double> contourPointPosition = std::make_pair(posi_x,posi_y);
            rawContour.push_back(contourPointPosition);
                if(debugMode==1){std::cout<<"for point "<<rawContour.size()<<" its posi is "<<contourPointPosition.first<<","<<contourPointPosition.second<<std::endl;};
        }
                if(debugMode==1){std::cout<<"File "<<inputFile<<" has "<<rawContour.size()<<" contour points"<<std::endl;};
        return rawContour;
    }

    void outputContour(const std::vector<std::pair<double,double>>& contour, const std::string& outputFileName){
        std::ofstream fout(outputFileName);
        fout<<"x y"<<std::endl;
        for(int i=0;i<contour.size();i++){
            fout<<contour[i].first<<" "<<contour[i].second<<std::endl;
        }
    }

    void vd_fout(const std::vector<double>& vd, const std::string& fileName){
        if(vd.size()==0){
            std::cerr<<"No elements in vd (vd_fout for"<<fileName<<")"<<std::endl;
            exit(1);
        }
        std::ofstream fout(fileName);
        for(int i=0;i<vd.size();i++){
            fout<<i<<" "<<vd[i]<<std::endl;
        }
        fout.close();
    }

    std::vector<std::string> getFilesNameInDirectory(const std::string& directory_path){
        std::vector<std::string> file_names;
        try{
            for(const auto& entry:std::filesystem::directory_iterator(directory_path)){
                if(entry.is_regular_file()){
                    file_names.push_back(entry.path().stem().string());
                }
            }
        }
        catch (const std::filesystem::filesystem_error& e){
            std::cerr<<"Filesystem error: "<<e.what()<<std::endl;
        } catch(const std::exception& e){
            std::cerr<<"Error: "<<e.what()<<std::endl;
        }
        return file_names;
    }

}