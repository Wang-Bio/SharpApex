#ifndef LINE2DV_H
#define LINE2DV_H

#include <vector>
#include <iostream>
#include "class2dv.h"

class Vertex;
class Organ;
class Line;
class Cell;

class Line{
 public:
    int vi[2];
    int li;
    std::vector<int> ci;
    Vertex d1;
    Vertex d2;
    bool IsOutermost;

    double length;
    double edgeForce;
    double frc_edge;
    int ordered_array;

    double slope, intercept;
    double A,B,C; //general form of a line: Ax+By=C

    double calc_length(Organ);
    double calc_length(void);
    void calc_slope_intercept(Organ);
    void calc_slope_intercept(void);
    void calc_general_ABC(Organ);
    double distance_from_point(Vertex);
    void set_endpoints(Vertex,Vertex);
    void set_endpoints(Vertex*,Vertex*);
    bool same_segment(Line);
    void print_slope_intercept(void);

    
    //~Line();
};



#endif