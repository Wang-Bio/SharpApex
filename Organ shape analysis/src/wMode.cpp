#include "../include/wMode.h"

void wMode::readFromFile(std::string modeFile){
        bool debugMode = 1;
    std::ifstream fin(modeFile,std::ios::in);
    std::string header;
    fin>>header;
    fin>>major;
    fin>>header;
    fin>>minor;
        if(debugMode==1){std::cout<<"Major mode is "<<major<<"; "<<"Minor mode is "<<minor<<std::endl;};
    fin.close();
}

void wMode::recordToFile(std::string modeFileOutput){
    std::ofstream fout(modeFileOutput);
    fout<<"Major_mode"<<std::endl;
    fout<<major<<std::endl;
    fout<<"Minor_mode"<<std::endl;
    fout<<minor<<std::endl;
    fout.close();
}