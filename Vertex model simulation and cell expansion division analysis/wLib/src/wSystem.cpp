#include "../include/wSystem.h"

wTimer global_timer;
wMode global_mode;

std::tm* wTimer::start(){
    start_time = std::chrono::high_resolution_clock::now();
    std::time_t current_time_t = std::chrono::system_clock::to_time_t(start_time);
    std::tm* local_time = std::localtime(&current_time_t);
    std::cout << "Program started at " << 1900 + local_time->tm_year 
    << "y " << 1 + local_time->tm_mon 
    << "m " << local_time->tm_mday 
    << "d " << local_time->tm_hour 
    << ": " << local_time->tm_min 
    << ": " << local_time->tm_sec << std::endl;
    running = true;

    std::ofstream fout(output_file);
    fout<< "Program started at " << 1900 + local_time->tm_year 
    << "y " << 1 + local_time->tm_mon 
    << "m " << local_time->tm_mday 
    << "d " << local_time->tm_hour 
    << ": " << local_time->tm_min 
    << ": " << local_time->tm_sec << std::endl;

    return local_time;
}

std::tm* wTimer::end(){
    if(running){
        end_time = std::chrono::high_resolution_clock::now();
        std::time_t current_time_t = std::chrono::system_clock::to_time_t(end_time);
        std::tm* local_time = std::localtime(&current_time_t);
        std::cout << "Program ended at " << 1900 + local_time->tm_year 
        << "y " << 1 + local_time->tm_mon 
        << "m " << local_time->tm_mday 
        << "d " << local_time->tm_hour 
        << ": " << local_time->tm_min 
        << ": " << local_time->tm_sec << std::endl;

        std::ofstream fout(output_file,std::ios::app);
        fout << "Program ended at " << 1900 + local_time->tm_year 
        << "y " << 1 + local_time->tm_mon 
        << "m " << local_time->tm_mday 
        << "d " << local_time->tm_hour 
        << ": " << local_time->tm_min 
        << ": " << local_time->tm_sec << std::endl;
        
        running = false;
        return local_time;
    }else{
        std::cerr<<"Timer is not running!"<<std::endl;
        return 0;
    }
}

void wTimer::totalRunningTime() const{
    if(!running){
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        int64_t total_milliseconds = duration.count();
        int hours = total_milliseconds/3600000;
        int minutes = (total_milliseconds%3600000)/60000;  
        int seconds = (total_milliseconds%60000)/1000;
        int miliseconds = total_milliseconds%1000;
        std::cout<<"This program lasted for "<<hours<< " hr,"<<minutes<<" min,"<<seconds<<" sec,"<<miliseconds<<" msecs."<<std::endl;
        std::ofstream fout(output_file,std::ios::app);
        fout<<"This program lasted for "<<hours<< " hr,"<<minutes<<" min,"<<seconds<<" sec,"<<miliseconds<<" msecs."<<std::endl;
    } else {
        std::cerr<< "Timer is still running!"<<std::endl;;
    }
}

void wMode::init(){
      std::ifstream fin(input_file,std::ios::in);
      if(!fin.is_open()){
         std::cerr<<"Error: missing "<<input_file<<" (the wMode file) !"<<std::endl;
         exit(1);
      }
      std::string tmp;
      fin>>tmp;
      fin>>major_mode;
      fin>>tmp;
      fin>>minor_mode;
      std::cout<<"Major mode : "<<major_mode<<" ; Minor mode: "<<minor_mode<<std::endl;
   }


namespace wSystem{
    void printMemoryUsage() {
        std::ifstream statusFile("/proc/self/status");
        std::string line;

        while (std::getline(statusFile, line)) {
            if (line.find("VmSize") != std::string::npos) {
                std::cout << line << std::endl;
            }
        }
    }

    std::vector<std::string> files_inside_folder(const std::string& folder){
        std::vector<std::string> files;
        for(const auto& entry : std::filesystem::directory_iterator(folder)){
            if(std::filesystem::is_regular_file(entry)){
                files.push_back(folder+entry.path().filename().string());
            }
        }
        return files;
    }
}