#ifndef WCONTOUR_H
#define WCONTOUR_H

#include <cmath>
#include <utility>
#include <vector>
#include <iostream>

#include "wFile.h"

namespace wContour{
    std::vector<std::pair<double,double>> normalization(const std::vector<std::pair<double,double>>&);
    std::vector<std::pair<double,double>> flipY(const std::vector<std::pair<double,double>>&); //imageJ get the y posi upside down, we need to flip it back to normal
    std::vector<std::pair<double,double>> parameterization(const std::vector<std::pair<double,double>>&, const int&);
    double calcPerimeter(const std::vector<std::pair<double,double>>&);

    std::vector<double> vec_curvature_circleFittingKasa_threeBoundaryPoints(const std::vector<std::pair<double,double>>&, const int&, const int&);
    double single_curvature_circleFittingKasa_threeBoundaryPoints(const std::pair<double,double>&, const std::pair<double,double>&, const std::pair<double,double>&);
    bool pointInPolygon_rayCasting(const std::pair<double,double>&, const std::vector<std::pair<double,double>>&);

    std::vector<double> vd_average(const std::vector<double>&, const int&);
    void vd_fout(const std::vector<double>&, const std::string&);
}   

#endif