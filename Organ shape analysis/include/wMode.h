#ifndef WMODE_H
#define WMODE_H

#include <string>
#include <iostream>
#include <fstream>

class wMode;

class wMode{
 private:
    std::string major;
    std::string minor;
 public:
    std::string get_major_mode(void){return major;};
    std::string get_minor_mode(void){return minor;};
    void set_major_mode(std::string mode_tmp){major=mode_tmp;};
    void set_minor_mode(std::string mode_tmp){minor=mode_tmp;};
    void readFromFile(std::string modeFileInput);
    void recordToFile(std::string modeFileOutput);
};

#endif