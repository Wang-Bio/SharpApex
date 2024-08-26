#ifndef _WSYSTEM_H
#define _WSYSTEM_H

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

class wTimer;
class wMode;

extern wTimer global_timer;
extern wMode global_mode;

class wTimer{
 private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time;
    bool running = false;
    std::string output_file = "simulation_timer.txt";
 public:
    std::tm* start();
    std::tm* end();
    void totalRunningTime() const;
};

class wMode{
 private:
   std::string major_mode;
   std::string minor_mode;
   std::string input_file = "../input/mode.txt";
 public:
   std::string get_major_mode(){
      return major_mode;
   }
   
   std::string get_minor_mode(){
      return minor_mode;
   }

   void init();

};

namespace wSystem{
   void printMemoryUsage(void);
   std::vector<std::string> files_inside_folder(const std::string&);
   
}

#endif