/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef GEO2DV_H
#define GEO2DV_H

#include "class2dv.h"
#include "vecInoue.h"
#include "parameter2dv.h"
#include "wangMath.h"
#include "force2dv.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <assert.h>
#include <tuple>

extern bool similarity_calculation_required;
extern vector<Vertex*> real_organ_contour_processed_for_similarity_index;
extern vector<tuple<int,double, double>> n_gons_lists; //int: n of n-gons, double: number of n-gons, double: the ratio of n-gons

namespace vertex_geo{
    double vertex_distance(Vertex*,Vertex*);
    double vertex_distance(Vertex,Vertex);
    bool vertex_relationship(Vertex*, Vertex*);
}

namespace line_geo{
    //lines
    bool line_vertex_relationship(Line,Vertex);
    pair<int,Vertex>  lines_relationship(const Line, const Line);
    _vec<double> line_cell_wall_intersection(double,double,Organ*,int);
    //segments
    double line_length(Organ*,int);
    Distance_point_line distance_line_segment_to_vertex(Organ*,int,Vertex);
    Distance_point_line distance_line_segment_to_vertex(Line,Vertex);
    
    bool line_segment_vertex_relationship(Line,Vertex);
    pair<int,Vertex> segments_relationship(Line,Line);

    vector<Line*> region_partition_lines_EdU(vector<Vertex*>);
    vector<double> angles_normalization_changing_axis(vector<Vertex*>,vector<Vertex*>);
}

namespace cell_geo{
    _vec<double> cell_center(Organ*, int);
    double cell_area(Organ*,int);
    double cell_perimeter(Organ*,int);
    double cell_regularity(Organ*,int);
    vector<int> cell_counterclock_line(Organ*,int);
}

namespace organ_geo{
    //1. geometric analysis of organ
    //1.1 basic 
    void organ_line_length(Organ*);
    double organ_cell_perimeter(Organ*);
    _vec<double> organ_center(Organ*);
    void epidermal_identity(Organ*);
    //1.2 area and perimeter
    double organ_area(Organ*);
    double organ_perimeter(Organ*);
    double organ_regularity(Organ*);
    double organ_circularity(Organ*);
    double organ_overlap(Organ*);
    //1.3 length width
    double organ_length(Organ*);
    double organ_width(Organ*);
    double organ_leaf_index(Organ*);
    _vec<double> organ_cell_arrangement(Organ*);

    double elliptical_index(Organ*);
    
    void organ_max_min_x_y(Organ*);

    //double organ_minimum_y_v(Organ*);
    //double organ_minimum_x_v(Organ*);
    //double organ_maximum_y(Organ*);
    //double organ_minimum_y(Organ*);

    //2. boundary analysis of organ
    //2.1 ordered boundary points generation
    vector<Vertex> organ_ordered_boundary_points_finding(Organ*,int);
    vector<Vertex*> organ_ordered_boundary_points_finding_pointer(Organ*,int);
    Ordered_boundary organ_ordered_anticlockwise_boundary(Organ*);
    vector<Vertex*> organ_boundary_points_along_polygon(Organ*,vector<Vertex*>,int);
    vector<Vertex> organ_boundary_points_along_polygon(Organ*,vector<Vertex>,int);
    vector<Vertex> organ_boundary_points_euclidean(Organ*,double);
    //2.2 curvature analysis
    vector<double> curvature_circle_fitting_kasa_three_boundary_points(Organ*,vector<Vertex*>,int,int);
    vector<double> curvature_circle_fitting_kasa_along_polygon(Organ*,int,int,int);
    vector<tuple<double,double,double>> organ_boundary_points_euclidean_point_number(Organ*,int);
    vector<double> curvature_arc_length_parameterization(Organ*,vector<Vertex*>,double);
    vector<double> curvature_graph_of_a_function(Organ*,vector<Vertex*>);

    double from_curvature_to_tip_position(Organ*,int,int,int);
    //3. complexity analysis of organ
    double complexity_index_information_theory(Organ*,vector<Vertex*>);
    pair<double,double> distance_entropy_error_information_theory(Organ*,int,vector<Vertex*>);
    //pair<double,double> local_angle_entropy_error(Organ*);
    double geometric_entropy(Organ*);

    //4. energy analysis of organ
    double organ_potential_energy(Organ*);

    //0 analysis tool
    bool point_in_polygon_ray_casting(Vertex*,Organ*);
    //bool point_in_polygon_ray_casting(Vertex,Organ*);
    vector<_vec<double>> line_polygon_intersection(double,double,Organ*);
    void organ_vertex_counterclockwise_sort(Organ*);

    //5. minor analysis
    vector<tuple<int,double,double>> n_gons_analysis(Organ*);
}

namespace geo{
    void calcGeometrics(Organ*);
    void calcGeometrics_cout(Organ*);

    void batch(Batch *, Organ *);

    void basic_VTK_analysis(Organ *);

    double distP(double x1, double y1, double x2, double y2);
    bool relaP(double x1, double y1, double x2, double y2);
    int relaL(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    pair<double, double> intersectLine(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    int relaRS(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
    //this ray is parallel to x axis and goes to infinite large 
    int relaRS_1(double xi1, double yi1, double xj1, double yj1, double xj2, double yj2);
    int relaS(double xi1, double yi1, double xi2, double yi2,double xj1, double yj1,double xj2, double yj2);
}

namespace geo_basic{
    typedef std::pair<_vec<double>, _vec<double>> pLine;
    bool isIntersectSegmentLine(pLine&, pLine&);
    extern double EPS;
    _vec<double> crossPoint(geo_basic::pLine &, geo_basic::pLine &);
}

namespace geo_vv{
    double maximum_distance_between_points_brute_force(vector<Vertex>,Organ*);
    double maxmum_distance_between_points_rotating_caliphers(vector<Vertex>,Organ*);
    vector<Vertex*> after_ImageJ_process(vector<Vertex*>);
    vector<Vertex*> normalization(vector<Vertex*>);
    vector<Vertex> normalization(vector<Vertex>);
    vector<Vertex> normalization_y_only(vector<Vertex>);
    bool equal_Euclidean_distance_except_last(vector<Vertex*>);
    vector<double> first_derivative(vector<Vertex*>);
    vector<double> second_derivative(vector<Vertex*>);
    vector<double> second_derivative_fourth_order(vector<Vertex*>);

    double area_vv_boundary(vector<Vertex*>);
    double perimeter_vv_boundary(vector<Vertex*>);
    double perimeter_vv_boundary(vector<Vertex>);
    vector<Vertex*> organ_boundary_points_along_polygon(vector<Vertex*>,int);
    vector<Vertex> organ_boundary_points_along_polygon(vector<Vertex>,int);
    bool point_in_polygon_ray_casting(Vertex*,vector<Vertex*>);
    bool point_in_polygon_ray_casting(Vertex,vector<Vertex>);
    vector<double> curvature_circle_fitting_kasa_three_boundary_points(vector<Vertex*>,int,int);
    vector<double> curvature_circle_fitting_kasa_three_boundary_points(vector<Vertex>,int,int);

    double vd_minimum(vector<double>);
    double vd_maximum(vector<double>);
    double accumulated_negative(vector<double>);
    vector<double> vd_averaged(vector<double>,int);

    vector<Vertex*> sort_vector_vertex_ascend(vector<Vertex*>);
    vector<Vertex*> sort_vector_vertex_descend(vector<Vertex*>);
    vector<Vertex> sort_vector_vertex_ascend(vector<Vertex>);
    vector<Vertex> sort_vector_vertex_descend(vector<Vertex>);
    bool comp_descend(double,double);
    double linear_fitting(Vertex*,Vertex*,double);
    double linear_fitting(Vertex,Vertex,double);
    vector<Vertex*> vv_x_swap(vector<Vertex*>);
    vector<Vertex> vv_x_swap(vector<Vertex>);
    vector<Vertex*> vector_vertex_sampling(vector<Vertex*> vv,double sampling_distance);
    vector<Vertex> vector_vertex_sampling(vector<Vertex> vv,double sampling_distance);
}

namespace boundary_geo{
    double surface_vertex_curvature_kasa(vector<Vertex*>);
    double surface_vertex_curvature_kasa(vector<Vertex>);
    double similarity_Index_1(vector<Vertex*>,vector<Vertex*>,double);
    double similarity_Index_1(vector<Vertex>,vector<Vertex*>,double);
    double similarity_Index_1(vector<Vertex>,vector<Vertex>,double);
    double similarity_Index_2(string,string,double);
    vector<Vertex*> read_and_process_real_organ_contour_imagej(string);
    double similarity_cal_during_simulation(Organ*,vector<Vertex*>);
    double similarity_for_contour_contour(std::string,std::string,double); //calculating the similarity of two .txt files with contour information 

    
}

namespace circle_geo{
    Intersection_relationship intersection_line_circle(Organ*,int,Circle);
    Intersection_relationship intersection_line_segment_circle(Organ*,int,Circle);
    Intersection_relationship intersection_line_segment_circle(Line,Circle);
}

namespace geo_analysis{
    void geo_plot(Geometrics_analysis_class, std::string, std::string, std::string);    
}



#endif