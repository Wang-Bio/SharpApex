#ifndef CELL2DV_H
#define CELL2DV_H

#include <vector>
#include <iostream>
#include "class2dv.h"

class Vertex;
class Line;
class Cell;
class Organ;

class Cell{
 public:
    std::vector<int> li;
    std::vector<int> vi;

    int n_edges;

    bool IsEpidermal;

    double area;
    double areaForce;

    int cellDivisionCount;
    double cellTime;

    double axisTheta;
    _vec<double> axis;
    double surfaceVertex[2];
    _vec<double> center;

    double outermostLength;
    
    double perimeter;

    int layer=-1;
    

    

    double regularity;
    double area_R;

    double y_rank;
    double area_rank;
    double frequency_modifier=1;

    int tag=0;

    double Gaussian_modifier=1;
    double area_modifier=1;

    double S_std;

    //~Cell();

};


#endif