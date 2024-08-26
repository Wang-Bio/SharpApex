/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _IO2DV_H
#define _IO2DV_H

#include "vecInoue.h"
#include "parameter2dv.h"
#include "class2dv.h"
#include "geo2dv.h"
#include "Initial2dv.h"

#include <fstream>
#include <iostream>
#include <dirent.h>
#include <unistd.h>
#include <cstring>
#include <ctime>
#include <tuple>
#include <cstdio>
#include <opencv2/opencv.hpp>

using namespace std;

namespace readV{
    void read_organ_txt(Organ*,int);
    
    void oneVTK(Organ*,string,string);
    void oneVTK(Organ*,string,string,int);
    void oneVTK(Organ*,std::pair<std::string,std::string>);
    void oneCell(Organ*,string);
    void oneLine(Organ*,string);
    void oneCellLine(Organ*);

    //void allVTK(Organ*,string);

    void final_VTK(Organ*,string);
    vector<Vertex*> xy_txt_to_vertexXY(string);
    vector<double> read_vd(string);
    vector<Vertex*> read_vv(string);
    vector<Vertex> read_vv_(string);

    Geometrics_analysis_class read_geo_output_txt(string);
}

namespace output{
    void VTK(Organ*);
    void VTK_debug(Organ*,std::string);
    
    void geo(Organ*);
    void division(Organ*);
    void batch(Batch*);
    void batch_final_analysis(vector<double>,vector<double>,vector<double>,string);
    void simulation_log(vector<time_t>, vector<time_t>);
    void geo_initial(void);

    void simulated_contour(Organ*, std::string);
    void simulated_contour(Organ*, std::string,std::string);
    void simulated_contour(Organ*);
}

namespace file_process{

int getFileNum(const string&);

}

namespace cout_fout_debug{
    void cout_vector_vertex(vector<Vertex*>);
    void cout_vector_vertex(vector<Vertex>);
    void fout_vector_vertex(vector<Vertex*>,string);
    void fout_vector_vertex(vector<Vertex>,string);
    void cout_vector_double(vector<double>);
    void cout_vector_int(vector<int>);
    void cout_vector_string(vector<string>);
    void cout_vector_line(vector<Line*>);

    void cout_vector_x_y(vector<double>,vector<double>);

    void fout_vector_double(vector<double>,string);
    void fout_vector_int_double(vector<int>,vector<double>,string);
    void fout_vector_double(vector<double>,vector<double>,vector<double>,string,string,string,string);
    void fout_vector_int(vector<int>,string);
    void fout_vector_x_y(vector<double>,vector<double>,string);

    void fout_pair_double(vector<pair<double,double>>,string);
    void fout_tuple_double(vector<tuple<double,double,double>>,string);
    void fout_tuple_double(vector<tuple<int,double,double>>,string);
    void fout_tuple_double(vector<tuple<int,int,double>>,string);


    void fout_heatmap_data(std::vector<int>,std::vector<int>,std::vector<double>,std::string);
    void fout_heatmap_data(std::vector<double>,std::vector<double>,std::vector<double>,std::string);
    std::vector<std::pair<double,double>> fin_vec_pair_double(std::string);
}



namespace gnu_plot{
    void dotplot(vector<double>,int,string);
    void organ_contour_plot(vector<Vertex*>,int,string);
    void organ_contour_plot(vector<Vertex*>);
    void curvature_plot(vector<double>,string);
    void panel_plot(vector<string>,int,int,string);
    void panel_plot_triangle(vector<string>,int,int,string);
    void png_to_avi(string);
    void gaussian_plot(double,double,string);

    void histogram(std::vector<double>,std::string);

    void multiple_time_lapse_contour_txt_plot(std::vector<std::string>, std::string, int);
    void organ_contour_plot(std::string,std::string,int);
    void organ_width_plot(std::string,std::string,int);

    void heatmap(std::string, std::string, int); 
    void heatmap(std::vector<int>, std::vector<int>, std::vector<double>, std::string); //not sure
    //void heatmap(std::vector<int>, std::vector<int>, std::vector<double>, std::string);
}

#endif