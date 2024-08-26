/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _CLASS2DV_H
#define _CLASS2DV_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "vecInoue.h"
#include "Vertex2dv.h"


constexpr int DEGREE_ACCURACY =2;
constexpr double EPS_geo = 1e-5;
using namespace std;
//list of class
class Organ;
class Vertex;
class Line;
class Cell;
class DivisionRecord;
class Geometrics;
class VertexCheck;
class CellCheck;
class Circle;
class Batch;
class parameterList;
class Intersection_relationship;
class Distance_point_line;
class Ordered_boundary;
class Quadratic_Solution;

//Quadratic_Solution Quadratic_Equation_Solve_(double,double,double);

class parameterList{
 public:
    std::vector<double> parameter_1;
    std::vector<double> parameter_2;
};

class Batch{
 public:
    std::vector<double> organ_perimeter;
    std::vector<double> organ_area;
    std::vector<double> averaged_perimeter;
    std::vector<double> averaged_inner_area;
    std::vector<double> averaged_epi_area;
    std::vector<int> N_in;
    std::vector<int> N_epi;
    std::vector<double> overlap_area;
    std::vector<double> real_area;
    std::vector<double> overlap_index;
    std::vector<double> averaged_regularity_in;
    std::vector<double> averaged_regularity_epi;
};

class Geometrics_analysis_class{
 public:
    std::vector<std::string> variable_name;
    std::vector<std::vector<double>> value;

    int findVariableIndex(const std::string &header){
        auto it = std::find(variable_name.begin(), variable_name.end(), header);
        if(it != variable_name.end()){
            return std::distance(variable_name.begin(), it);
        }
        return -1;
    }

    std::vector<double> getColumnData(std::string header){
        int idx = findVariableIndex(header);
        if(idx == -1){
            std::cerr<<"Error: Variable name not found!"<<std::endl;
            return {}; //return an empty vector
        }

        std::vector<double> column;
        for(const auto& row : value){
            if(idx<row.size()){
                column.push_back(row[idx]);
            }
        }
        return column;
    }
};


class DivisionRecord{
 public:
    int time;
    int cidx;
    bool IsEpidermal;
    double axisTheta; 

    double center_x;
    double center_y;

    int division_count;
    double av_in_frequency_modifier;
    double av_epi_frequency_modifier;
};

class Circle{
 public:
    _vec<double> center;
    double radius;
};

class Quadratic_Solution{
 public:
    double delta;
    double x1;
    double x2;
};

class Intersection_relationship{
 public:
    int Relationship;
    vector<Vertex> cross_points;
};

class Distance_point_line{
 public:
    double distance;
    Vertex Closest_Point;
    double t;
};

class Ordered_boundary{
 public:
    vector<int> li;
    vector<Vertex> vi;
};

/*
Quadratic_Solution Quadratic_Equation_Solve_(double A_tmp,double B_tmp,double C_tmp)
    {
        //quadratic equation: 𝐴𝑥^2+𝐵𝑥+𝐶=0
        //The solution is 𝑥_1,2=(−𝐵±√(𝐵^2−4𝐴𝐶))/2𝐴, 𝑖𝑓 √(𝐵^2−4𝐴𝐶)≥0
        Quadratic_Solution Result_Quadratic;
        Result_Quadratic.delta = B_tmp*B_tmp - 4*A_tmp*C_tmp;

        if(Result_Quadratic.delta>0){
            Result_Quadratic.x1 = (-B_tmp+sqrt(Result_Quadratic.delta))/(2*A_tmp);
            Result_Quadratic.x2 = (-B_tmp-sqrt(Result_Quadratic.delta))/(2*A_tmp);
        }
        else if(Result_Quadratic.delta==0){
            Result_Quadratic.x1 = -B_tmp/(2*A_tmp);
        }
        else{

        }

        return Result_Quadratic;
    }
*/

#endif
