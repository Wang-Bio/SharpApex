#ifndef ORGAN2DV_H
#define ORGAN2DV_H

#include <vector>
#include "vecInoue.h"
#include "class2dv.h"

class Vertex;
class Line;
class DivisionRecord;
class Cell;
class Geometrics;
class ForceCheck;

class Organ{
 public:
    int step;
    std::vector<Vertex *> p_v;
    std::vector<Line *>p_l;
    std::vector<Cell *>p_c;
    std::vector<DivisionRecord *> d_r;

    std::vector<ForceCheck* >forceCheck;
    double initial_cell_number;

    //geometrics analysis
    //1. epidermal identity
    int N_inner_cell;
    int N_epi_cell;
    std::vector<int> surface_vertex;
    std::vector<int> surface_line;
    //2. organ center
    _vec<double> center;
    //3. organ cell layer
    int cell_layer_number;
    //4. organ area
    double area; //the summed area for all cells
    double epiArea; //the summed area for all epidermal cells
    double inArea;  //the summed area for all inner cells
    double epiArea_averaged;
    double inArea_averaged;
    double area_averaged;
    //5. organ perimeter
    double perimeter;
    double perimeter_averaged;
    //6. organ circularity
    double circularity; 
    //7. organ regularity
    double regularity_averaged;
    double regularity_in_av;
    double regularity_epi_av;
    //8. organ length, width
    double organ_length;
    double organ_width;
    std::pair<_vec<double>,_vec<double>> length_vertex;
    std::pair<_vec<double>,_vec<double>> width_vertex;
    double length_axis_slope;
    double length_axis_intercept;
    double width_axis_slope;
    double width_axis_intercept;
    double leaf_index;
    double length_width_ratio;
    //9. cell arrangement
    double cell_arrangement_x;
    double cell_arrangement_y;
    double cell_arrangement_ratio;
    //10. cell elliptical index
    double elliptical_index;
    //11. organ overlap
    std::vector<Vertex*> intersection;
    std::vector<Vertex*> intersection_real;
    std::vector<Vertex*> surface_pv_real;
    std::vector<Line*> surface_pl_real;
    double realArea;
    double overlapArea;
    double overlap_index;
    //12. similarity index
    double similarity_index;
    //13. curvature
    double maximum_curvature;
    double minimum_curvature;
    double accumulated_negative_curvature;

    //analysis of area control division frequency ratio between epidermal cells and inner cells
    double av_in_division_frequency_modifier;
    double av_epi_division_frequency_modifier;
    double F_ratio; //divsion frequency ratio between epidermal cells and inner cells
    

    //more geometrics
    double organ_potential_energy;
    double y_min_vertex;
    double y_max_vertex;
    double x_min_vertex;
    double x_max_vertex;
    double y_min_cell;
    double y_max_cell;
    double x_min_cell;
    double x_max_cell;

    double tip_length_min_c_position_absolute;
    double tip_length_min_c_position_relative;

    ~Organ();

    //debugs
    void print_basics(void){
        std::cout<<"Cell number: "<<p_c.size()<<"; line number: "<<p_l.size()<<"; vertex number: "<<p_v.size()<<std::endl;
    }
};



#endif