#ifndef WANGSYSTEM_H
#define WANGSYSTEM_H
#include <stdio.h>
#include <math.h>

#include <ctime>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>

using namespace std;
 
extern std::vector<time_t> initial_time_v, terminal_time_v; 

namespace wangSystem{
    time_t initial_time(void);
    time_t terminal_time(void);
    void printMemoryUsage(void);
    std::vector<std::string> files_inside_folder(const std::string&);
}

bool areAllValuesGreaterThan(std::vector<double>&,int,double);


#endif