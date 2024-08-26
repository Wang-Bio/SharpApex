#ifndef WBMP_H
#define WBMP_H

#include <vector>
#include <string>

#include <iostream>
#include <fstream>
#include <cmath>

#include "wPair.h"
#include "wFile.h"

namespace wBMP{
    int rangeTics(int);
    void drawSmoothCircle(std::vector<uint8_t>&, int, int, int, int, int, uint8_t, uint8_t, uint8_t);
    void drawLine(std::vector<uint8_t>&, int, int, int, int, int, int, uint8_t, uint8_t, uint8_t);
    std::vector<int> generateTicks(int,int);
    void drawContour(std::vector<std::pair<double,double>>,std::string);
    void drawCurvature(std::vector<double>,std::string);
}

#endif