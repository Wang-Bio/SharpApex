#ifndef VERTEX2DV_H
#define VERTEX2DV_H

#include "vecInoue.h"
#include "class2dv.h"

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

class Vertex;
class Line;
class Cell;
class Organ;

class Vertex{
 public:
    std::vector<int> li; 
    std::vector<int> ci;
    _vec<double> loc; //Cartesian coordinates
    _vec<double> frc_area;
    _vec<double> frc_edge;
    
    double r,theta; //Polar coordinates
    unsigned int occurrenceInCell;
    bool IsSurface;

    std::vector<int> layer;
    int vi;
    std::vector<int> li_array;
    int vi_array;
    //destructor
   
    void print_Cartesian(void);
    void print_Polar(void);
    double distance_from_vertex(Vertex);
    void Cartesian_to_Polar(void);
    void Polar_to_Cartesian(void);
    bool collinear_points(Vertex,Vertex);
    bool same_vertex(Vertex);

};

namespace vVertex{

    std::vector<Vertex> read_from_txt(const std::string &,int,int);
    void save_txt(const std::string &, const std::vector<Vertex> &);
    
    void print(const std::vector<Vertex> &);
    std::vector<Vertex> vertical_reflection(const std::vector<Vertex> &);
    std::vector<Vertex> normalization(const std::vector<Vertex> &);
    std::vector<Vertex> sample(const std::vector<Vertex> &);
    std::vector<Vertex> vec_averaged(const std::vector<std::vector<Vertex>> &);
    int index_for_y(const std::vector<Vertex> &,double);
    void dot_plot(const std::vector<Vertex> &);
}

#endif