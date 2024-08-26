#ifndef ANALYSIS2DV_H
#define ANALYSIS2DV_H

#include "class2dv.h"
#include "vecInoue.h"
#include "parameter2dv.h"
#include "geo2dv.h"
#include "Organ2dv.h"

class analysis_2D;

class analysis_2D{
 private: 
    std::vector<std::pair<double,double>> param_2D;
    std::string param_name_1;
    std::string param_name_2;

    int repeat_time;
    std::vector<std::string> file_directory;
    
    //geometrics analysis
    //1. curvature analysis
    std::vector<double> min_curvature_all;
    std::vector<double> min_curvature_av;
    std::vector<double> max_curvature_all;
    std::vector<double> max_curvature_av;
    std::vector<double> range_curvature_all;
    std::vector<double> range_curvature_av;
    std::vector<double> accumulated_negative_curvature_all;
    std::vector<double> accumulated_negative_curvature_av;
    //2. similarity analysis
    std::vector<double> similarity_index_all;
    std::vector<double> similarity_index_av;
 
 public:
    //defualt constructor
    analysis_2D(){}

    //add parameters data
    void add_param_name(std::string, std::string);
    void add_param_single(double,double);
    void add_param_rect(std::vector<double>,std::vector<double>);
    void add_param_rect_incremental(double,double,double,double,double,double);
    void add_param_triang_right_up(std::vector<double>,std::vector<double>);
    void add_param_triang_right_up_incremental(double,double,double,double,double,double);
    
    //file_directory
    void add_all_file_directory(std::string);
    void add_all_file_directory_3(std::string);

    std::string get_file_directory(int);
    std::pair<std::string, std::string> VTK_directory(int,int);
    std::pair<std::string, std::string> VTK_directory(int);

    std::string get_parameter_s(int);

    //basic_properties
    void set_repeat_time(int);
    int get_repeat_time(void);
    int size(void);
    std::pair<double,double> get_parameter(int);
    std::pair<std::string,std::string> get_param_name(void);

    //output
    void print(void);

    //geometrics_recording
    void geometrics_recording(Organ*,int,int);
    void geometrics_output(std::string);

    void add_curvature_analysis(Organ*,int,int);
    void add_similarity_analysis(Organ*,int,int);
};

#endif