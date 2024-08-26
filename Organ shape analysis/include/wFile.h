#ifndef WFILE_H
#define WFILE_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <string>
#include <vector>
#include <filesystem>

extern std::string inputFileDirectory;
extern std::string outputFileDirectory;

namespace wFile{
    std::vector<std::pair<double,double>> readVecPairDouble(const std::string&);

    std::vector<std::pair<double,double>> readContour_posi_only(const std::string&);
    std::vector<std::pair<double,double>> readContour_raw(const std::string&);

    std::vector<std::string> getFilesNameInDirectory(const std::string&);

    void outputContour(const std::vector<std::pair<double,double>>&, const std::string&);
    void vd_fout(const std::vector<double>&, const std::string&);

    template <typename T>
    void vec_cout(const std::vector<T>& vec){
        for(int i=0;i<vec.size();i++){
            std::cout<<i<<" "<<vec[i]<<std::endl;
        }
    }

}

#endif