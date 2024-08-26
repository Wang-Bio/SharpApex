#include "../include/wangSystem.h"

std::vector<time_t> initial_time_v, terminal_time_v; 

namespace wangSystem{

    time_t initial_time(void){
        time_t curr_time;
        tm * curr_tm;
        char initial_time_string [100];
        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(initial_time_string, 100, "Simulation starts at %T, %Y %B %d", curr_tm);
        cout << initial_time_string << endl;
        initial_time_v.push_back(curr_time);
        return curr_time;
        //reference: https://www.programiz.com/cpp-programming/library-function/ctime/strftime
    }

    time_t terminal_time(void){
        time_t curr_time;
        tm * curr_tm;
        char terminal_time_string [100];
        time(&curr_time);
        curr_tm = localtime(&curr_time);
        strftime(terminal_time_string, 100, "Simulation terminates at %T, %Y %B %d", curr_tm);
        cout << terminal_time_string << endl;
        terminal_time_v.push_back(curr_time);
        return curr_time;
        //reference: https://www.programiz.com/cpp-programming/library-function/ctime/strftime
    }

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

bool areAllValuesGreaterThan(std::vector<double>& vec, int startIndex, double threshold) {
    // Check if startIndex is within the range of the vector
    if (startIndex < 0 || startIndex >= vec.size()) {
        return false; // or throw an exception
    }

    return std::all_of(vec.begin() + startIndex, vec.end(), [threshold](double value) {
        return value > threshold;
    });
}