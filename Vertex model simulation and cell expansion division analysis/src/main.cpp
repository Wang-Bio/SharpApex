/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <vector>
#include <random>
#include <time.h>
#include <cstring>
#include <assert.h>
#include <ctime>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <opencv2/opencv.hpp>
#include <filesystem>
#include <experimental/filesystem>

#include "../include/class2dv.h"
#include "../include/division.h"
#include "../include/vecInoue.h"
#include "../include/parameter2dv.h"
#include "../include/force2dv.h"
#include "../include/geo2dv.h"
#include "../include/IO2dv.h"
#include "../include/wangSystem.h"
#include "../include/wangMath.h"

#include "../include/vtk2dv.h"
#include "../include/gnuplot.h"
#include "../include/image2dv.h"

#include "../include/Organ2dv.h"
#include "../include/Cell2dv.h"
#include "../include/Line2dv.h"
#include "../include/Vertex2dv.h"
#include "../include/Initial2dv.h"
#include "../include/Analysis2dv.h"

#include "../wLib/include/wUtility.h"
#include "../wLib/include/wSystem.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkPNGWriter.h>
#include <vtkProperty.h>
#include <cstring>
#include <vtkRenderWindowInteractor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkLight.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkConeSource.h>
#include <vtkCellData.h>
#include <vtkTimerLog.h>


//Cell division frequency control:
//1. No control "no"; 2. "area"; 3. "balance"; 4. "Gaussian_area"; 5. "Gaussian balance"; 6. "balance_BF_3";
//7. "balance_BF_1"; 8. "temporal_Gaussian_linear_mu"; 9. "temporal_Gaussian_linear_sigma"
//10. "random_pick"; 11. "Gaussian_xy"; 12. "uniform_to_Gaussian"; 13. "biregion_frequency_position"; 14. biregion_frequency_identity
//15. "arrest_front" //16. "arrest_front_biregion"


//Cell division angles control:
//1. "random"; 2."anticlinal" (only for epidermal cells); 3. "periclinal" (only for epidermal cells); 4. "constant_0"; 5. "constant_90"
//6. "mochizuki_bias"; 7. "mochizuki_bias_asymmetrical"; 8. "yin_distribution"; 9."mochizuki_bias_apical_basal"; 10. "temporal_angle_bias"
//11. "Gaussian_bias"; 12. "biregion_angles_position"; 13. "Gaussian_bias_continuous_Gaussian"; 14. biregion_angles_identity
//15. "temporal_biregion_angles" 16. "biregion_angles_gradual"

//Cell division frequency and angles combination control: 
//1. "biregion_frequency_angles_position"

//Mechanics mode 
//1. "simple"; 2. L_std; 3. "sigma_O_spatial"; 4. "spatial_S_std"

using namespace std;

//string major_mode = "test"; 
//Major Mode: analysis, simualtion, test, or preparation

//string minor_mode = "single"; 

//minor mode for simulation
//single: only for a single simulation; repeat: repeat the simulation with same parameter for several times; batch: do simulation with defined chaning parameters

//minor mode for analysis


//initial: initialization txt file; single : single time frame; final: the final time frame of a single simulation; batch_final: the final time frame of batch simulations
//organ outline: txt_outline_analysis


int main(){

parameter::read_mode(modeFile);


//recording time
vector<time_t> initial_time, terminal_time;
time_t initial_time_first = wangSystem::initial_time();
initial_time.push_back (initial_time_first);

cout<<"******************| Mode Selection| ********************************"<<endl;
cout<<"Major mode: "<<major_mode<<endl;
if(major_mode=="simulation"||major_mode=="analysis"||major_mode=="experiment"||major_mode=="plot")
cout<<"Minor mode: "<<minor_mode<<endl;
cout<<"********************************************************************"<<endl;


if(major_mode=="simulation"){

if(minor_mode=="single"){

//initialization
//parameter::read_and_record();
parameter::read(parameterInputFile);
parameter::record();

Organ *p_g = new Organ;
readV::read_organ_txt(p_g,0);
initialization::organ(p_g);

//simulation
cout<<"*********************| Simulation Start |***************************"<<endl;
for(int step=0; step<step_end; step++){
    
    p_g->step=step;
    force::calcForceMotion(p_g);

    //check cell division
    if(step%T_division_check==0){
        cout<<"********************* At time step "<<step<<"***************************"<<endl;
        //cout<<"Start geometrics calculation ---- ";
        geo::calcGeometrics(p_g);
        //cout<<"Output ---- ";
        output::geo(p_g);
        //cout<<"Finished "<<endl;

        //cout<<"Start division event judgement ---- ";
        division::Global(p_g);
        //cout<<"Division events output ---- ";
        output::division(p_g);
        //cout<<"Finished "<<endl;
        //autoVTK::VTKLineCell(p_g);
        geo::calcGeometrics(p_g);
        cout<<"********************************************************************"<<endl;
    }

    //cell time add
    if(step%T_cell_time==0){
        division::cellTimeAdd(p_g,1);
    }

    //output VTK
    if(step%T_vtkoutput==0){
        output::VTK(p_g);            
        output::simulated_contour(p_g); 
        if(termination::checkEnd(p_g)==1){
            goto End_One_Simulation1;
        }
    }
}
End_One_Simulation1:;
    cout<<"**********************| Simulation End |****************************"<<endl;

}

else if(minor_mode=="repeat"){

int repeat_time = 3;
parameter::read(parameterInputFile);
parameter::record();
initialOrganFile= "../"+initialOrganFile;

for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){

char dir_file[100];
sprintf(dir_file,"%d/",repeat_i);
mkdir(dir_file,0777);
chdir(dir_file);

//initialization
Organ *p_g = new Organ;
readV::read_organ_txt(p_g,0);
initialization::organ(p_g);

//simulation
cout<<"*********************| Simulation Start |***************************"<<endl;
for(int step=0; step<step_end; step++){
    p_g->step=step;
    force::calcForceMotion(p_g);

    //check cell division
    if(step%T_division_check==0){
        cout<<"At time step "<<step<<endl;
            //cout<<"Calculate geometrics information ---- ";
            geo::calcGeometrics(p_g);
            //cout<<"Output ---- ";
            output::geo(p_g);
            //cout<<"Finished "<<endl;

            //cout<<"Start division events judgement ---- "<<endl;
            division::Global(p_g);
            //cout<<"Division events output ---- ";            
            output::division(p_g);
            //cout<<"Finished "<<endl;
            geo::calcGeometrics(p_g);
        cout<<"********************************************************************"<<endl;
    }

    //cell time add
    if(step%T_cell_time==0){
        division::cellTimeAdd(p_g,1);
    }

    //output VTK
    if(step%T_vtkoutput==0){
        output::VTK(p_g);    
        output::simulated_contour(p_g);
        if(termination::checkEnd(p_g)==1){
                    goto End_One_Simulation2;
                }
    }
}
End_One_Simulation2:;
    cout<<"********************************************************************"<<endl;
chdir("../");
}

}

else if(minor_mode=="batch"){
int repeat_time =3;
string parameter_1_str = "relative frequency";
string parameter_1_str_abbr = "rf"; //for creating filefolders
string parameter_2_str = "beta apical";
string parameter_2_str_abbr = "ba"; //for creating filefolders

parameterList* pL = new parameterList;
parameter::batchRead(pL);
//parameter::read("../"+parameterInputFile);
parameter::read(parameterInputFile);
initialOrganFile= "../"+initialOrganFile;

for(int pL_i=0;pL_i<pL->parameter_1.size();pL_i++){
    Batch *ba = new Batch;
    biregion_frequency_angles_position_relative_frequency = pL->parameter_1[pL_i];
    biregion_frequency_angles_position_apical_gaussian_beta = pL->parameter_2[pL_i];
    if(in_division_direction=="random"){
        cout<<"Current "<<division_control<<" control "<<parameter_1_str<<" "<<parameter_2_str<<": "<<pL->parameter_1[pL_i]<<" "<<pL->parameter_2[pL_i]<<endl;
    }
    else{
        cout<<"Current "<<in_division_direction<<" control "<<parameter_1_str<<" "<<parameter_2_str<<": "<<pL->parameter_1[pL_i]<<" "<<pL->parameter_2[pL_i]<<endl;
    }
    char dir_file[100];
    sprintf(dir_file,"%s_%.2f_%s_%.2f/",parameter_1_str_abbr.c_str(),pL->parameter_1[pL_i],parameter_2_str_abbr.c_str(),pL->parameter_2[pL_i]);
    mkdir(dir_file,0777);
    chdir(dir_file);

    parameter::record();

    time_t initial_time_tmp = wangSystem::initial_time();
    initial_time.push_back(initial_time_tmp);
for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
    char dir_file[100];
    sprintf(dir_file,"%d/",repeat_i);
    mkdir(dir_file,0777);
    chdir(dir_file);

    //initialization
    Organ *p_g = new Organ;
    readV::read_organ_txt(p_g,0);
    organ_geo::epidermal_identity(p_g);

    force::forceShapeInitiation(p_g,200000);
    division::cell_time_initialization(p_g);

    //simulation
    cout<<"*********************| Simulation Start |***************************"<<endl;
    for(int step=0; step<step_end; step++){
        p_g->step=step;
        force::calcForceMotion(p_g);

        //check cell division
        if(step%T_division_check==0){
            cout<<"At time step "<<step<<endl;
            //cout<<"Calculate geometrics information ---- ";
            geo::calcGeometrics(p_g);
            //cout<<"Output ---- ";
            output::geo(p_g);
            //cout<<"Finished "<<endl;

            //cout<<"Start division events judgement ---- "<<endl;
            division::Global(p_g);
            //cout<<"Division events output ----  ";            
            output::division(p_g);
            //cout<<"Finished "<<endl;
            geo::calcGeometrics(p_g);
            cout<<"********************************************************************"<<endl;
        }

        //cell time add
        if(step%T_cell_time==0){
            division::cellTimeAdd(p_g,1);
        }

        //output VTK
        if(step%T_vtkoutput==0){
            output::VTK(p_g);     
            output::simulated_contour(p_g);       
            if(termination::checkEnd(p_g)==1){
                        goto End_One_Simulation3;
                    }
        }
    }
    End_One_Simulation3:;
    cout<<"********************************************************************"<<endl;
    chdir("../");
    }
    time_t terminal_time_tmp = wangSystem::terminal_time();
    terminal_time.push_back(terminal_time_tmp);
    chdir("../");
}
}

else if(minor_mode=="continue"){
    Organ *p_g = new Organ;

    string FileDirectory = "../../basic_simulations_results/nankinhaze/";
    //string FileDirectory = "/mnt/d/c_file/2023/2023_02_23_biregion_angles_identity/y_0.60/y_0.60_ab_20.00/0/";
    //initialization
    //parameter::read_and_record();
    parameter::read(parameterInputFile);
    chdir(FileDirectory.c_str());
    parameter::record();
    readV::final_VTK(p_g, FileDirectory);
    initialization::organ_continue(p_g);
    p_g->step++;
    //simulation
    cout<<"*********************| Simulation Start |***************************"<<endl;
    for(int step=(p_g->step*T_vtkoutput); step<step_end; step++){
        
        p_g->step=step;
        force::calcForceMotion(p_g);

        //check cell division
        if(step%T_division_check==0){
            cout<<"At time step "<<step<<endl;
            //cout<<"Start geometrics calculation ---- ";
            geo::calcGeometrics(p_g);
            //cout<<"Output ---- ";
            output::geo(p_g);
            //cout<<"Finished "<<endl;

            //cout<<"Start division event judgement ---- ";
            division::Global(p_g);
            //cout<<"Division events output ---- ";
            output::division(p_g);
            //cout<<"Finished "<<endl;
            //autoVTK::VTKLineCell(p_g);
            cout<<"********************************************************************"<<endl;
        }

        //cell time add
        if(step%T_cell_time==0){
            division::cellTimeAdd(p_g,1);
        }

        //output VTK
        if(step%T_vtkoutput==0){
            output::VTK(p_g);        
            output::simulated_contour(p_g);    
            if(termination::checkEnd(p_g)==1){
                        goto End_One_Simulation4;
                    }
        }
    }
    End_One_Simulation4:;
        cout<<"********************************************************************"<<endl;


    }

    else{
        cout<<"Fatal error: no minor mode selected ! (major mode is simulation)"<<endl;
        exit(-1);
    }

}

else if(major_mode=="analysis"){

if(minor_mode=="single"){

Organ* p_g = new Organ;

string Cell_VTK_Filename ="/mnt/e/c_file/2023_09_13_biregion_angles_gradual/y1_0.7_y2_0.5/0/2dv_cell0018900000.vtk";
string Line_VTK_Filename ="/mnt/e/c_file/2023_09_13_biregion_angles_gradual/y1_0.7_y2_0.5/0/2dv_line0018900000.vtk";
//int boundary_points_number = 100;
//int fitting_distance=10;
//int NumAveraging = 9;

//************************Analyzing cell area distribution***************************************//

readV::oneVTK(p_g, Cell_VTK_Filename, Line_VTK_Filename);
std::cout<<"Cell number "<<p_g->p_c.size()<<std::endl;
geo::calcGeometrics(p_g);

int boundary_points_number = 100;
int fitting_distance = 10;
int NumAveraging = 9;
organ_geo::from_curvature_to_tip_position(p_g,boundary_points_number,fitting_distance,NumAveraging);

//cout_fout_debug::fout_vector_double(curvature_kasa,"TS_new_curvature.txt");
//gnu_plot::dotplot(curvature_kasa,1,"TS_new_curvature.png");


//std::cout<<"similarity_index "<<p_g->similarity_index<<std::endl;
//std::cout<<"Minimum curvature "<<p_g->minimum_curvature<<std::endl;
/*
double sampling_distance = 0.001;
double similarity_index;
vector<Vertex> simulated_contour;
for(int vi=0;vi<(int)p_g->p_v.size();vi++){
    if(p_g->p_v[vi]->IsSurface==1){
        simulated_contour.push_back(*p_g->p_v[vi]);
    }
}
vector<Vertex> simulated_contour_normalized_y_only = geo_vv::normalization_y_only(simulated_contour);
//vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
vector<Vertex> simulated_outline_sampled_y_only = geo_vv::vector_vertex_sampling(simulated_contour_normalized_y_only,sampling_distance);
//cout_fout_debug::cout_vector_vertex(simulated_outline_sampled_y_only);

vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
//vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);

std::vector<double> contour_widths;
for(int width_i=0;width_i<1/sampling_distance;width_i++){
    double contour_width = abs(simulated_outline_sampled_y_only[width_i].loc.x) + abs(simulated_outline_sampled_y_only[2/sampling_distance-width_i].loc.x);
    contour_widths.push_back(contour_width);
}
cout_fout_debug::fout_vector_vertex(simulated_outline_sampled, "contour_2000_cell.txt");
cout_fout_debug::fout_vector_double(contour_widths, "contour_absolute_width_2000_cell.txt");

vector<Vertex> basal_contour; //tip removed
vector<int> tip_removed_index;
    //the very tip should be directly removed
    tip_removed_index.push_back(1/sampling_distance);

double absolute_width_threshold = 25.0;
for(int i=(int)1/sampling_distance-1; i>0; i-=sampling_distance){
    if(contour_widths[i]<absolute_width_threshold){
        tip_removed_index.push_back(i);
        tip_removed_index.push_back(2/sampling_distance-i);
    }
    else{
        break;
    }
}

for(int i=0;i<simulated_outline_sampled.size();i++){
    bool is_true=1;
    for(int j=0;j<tip_removed_index.size();j++){
        if(i==tip_removed_index[j]){
            is_true=0;
        }
    }

    if(is_true==1){
        basal_contour.push_back(simulated_outline_sampled[i]);
    }
}
cout_fout_debug::cout_vector_int(tip_removed_index);
cout_fout_debug::fout_vector_vertex(basal_contour, "contour_basal_2000_cell.txt");
*/
//vector<Vertex> simulated_outline_swapped_sampled = geo_vv::vector_vertex_sampling(simulated_contour_swapped,sampling_distance);
/*
std::vector<double> cell_area_all;
double cell_area_av=0;
int cell_number_for_cell_area_av=0;
for(int ci=0;ci<(int)p_g->p_c.size();ci++){
    if(p_g->p_c[ci]->IsEpidermal==1){
        cell_area_all.push_back(p_g->p_c[ci]->area);
        cell_area_av += p_g->p_c[ci]->area;
        cell_number_for_cell_area_av++;
    }
}
cell_area_av /= cell_number_for_cell_area_av;
double cell_area_va =0;
for(int ci=0;ci<(int)p_g->p_c.size();ci++){
    if(p_g->p_c[ci]->IsEpidermal==1){
        cell_area_va += std::pow((p_g->p_c[ci]->area-cell_area_av),2);
    }
}
cell_area_va/=cell_number_for_cell_area_av;

std::cout<<"cell area averaged: "<<cell_area_av<<"; cell area variance: "<<cell_area_va<<std::endl;

gnu_plot::histogram(cell_area_all,"nankinhaze_shape_epi_cell_area_distribution.png");
*/
/*
geo::basic_VTK_analysis(p_g);

geo::calcGeometrics(p_g);
vector<Vertex> anticlockwise_surface_vertices = organ_geo::organ_ordered_anticlockwise_boundary(p_g).vi;
anticlockwise_surface_vertices = geo_vv::normalization(anticlockwise_surface_vertices);
//cout_fout_debug::fout_vector_vertex(anticlockwise_surface_vertices,"../analysis_example/petal_contour.txt");
vector<Vertex> boundary_points = organ_geo::organ_ordered_boundary_points_finding(p_g, boundary_points_number);
boundary_points = geo_vv::normalization(boundary_points);
//cout_fout_debug::fout_vector_vertex(boundary_points,"../analysis_example/petal_boundary_equal_distance.txt");
vector<double> curvature_kasa = organ_geo::curvature_circle_fitting_kasa_along_polygon(p_g,boundary_points_number,fitting_distance,NumAveraging);
cout_fout_debug::fout_vector_double(curvature_kasa,"../analysis_result/petal_curvature.txt");
gnu_plot::dotplot(curvature_kasa,1,"../analysis_result/petal_curvature_example.png");

cout<<"elliptical index: "<<organ_geo::elliptical_index(p_g)<<endl;
cout<<"circularity index: "<<p_g->circularity<<endl;
cout<<"Geometric entropy: "<<organ_geo::geometric_entropy(p_g)<<endl;
double minimum_curvature = geo_vv::vd_minimum(curvature_kasa);
double accumulated_curvature = geo_vv::accumulated_negative(curvature_kasa);
cout<<"minimum_curvature: "<<minimum_curvature<<"; accumulated_curvature: "<<accumulated_curvature<<endl;
*/
/*

//Ordered_boundary ordered_boundary_tmp = organ_geo::organ_ordered_anticlockwise_boundary(p_g);
//cout_fout_debug::fout_vector_vertex(ordered_boundary_tmp.vi,"ordered_boundary_tmp.txt");
//double boundary_point_distance = 1.0;
//vector<Vertex*> boundary_points = organ_geo::organ_boundary_points_euclidean(p_g,boundary_point_distance);
//cout<<"Does boundary points have the same distance ? : "<<geo_vv::equal_Euclidean_distance_except_last(boundary_points)<<endl;
//cout<<"Boundary point distance "<<boundary_point_distance<<"; number of boundary points: "<<boundary_points.size()<<endl;
//cout<<"The whole arc-length of the curve: "<<boundary_point_distance*boundary_points.size()<<endl;
//vector<double> curvature_boundary = organ_geo::curvature_arc_length_parameterization(p_g,boundary_points,boundary_point_distance);
//cout_fout_debug::fout_vector_vertex(boundary_points,"boundary_points_distance_1.0_higher_order.txt");
//cout_fout_debug::fout_vector_double(curvature_boundary,"boundary_curvature_distance_1.0_higher_order.txt");

//vector<double> curvature_boundary_graph_function = organ_geo::curvature_graph_of_a_function(p_g, boundary_points);
//cout_fout_debug::fout_vector_double(curvature_boundary_graph_function,"boundary_curvature_distance_1.0_graph_function.txt");

//organ_geo::organ_boundary_points_euclidean_point_number(p_g,100);

//organ_geo::organ_length(p_g);
//organ_geo::organ_width(p_g);
//organ_geo::organ_leaf_index(p_g);
//organ_geo::organ_cell_arrangement(p_g);
//cout<<"organ length: "<<p_g->organ_length<<"; organ width: "<<p_g->organ_width<<"; organ leaf index: "<<p_g->organ_length/p_g->organ_width<<endl;
//cout<<"organ cell arrangement x: "<<p_g->cell_arrangement_x<<" ; organ cell arrangement y: "<<p_g->cell_arrangement_y<<" ; organ cell arrangement ratio: "<<p_g->cell_arrangement_ratio<<endl;
//output::geo(p_g);
*/

}

else if(minor_mode=="single_time_lapse"){
    string FileDirectory = "/mnt/e/c_file/2023_09_18_arrest_front_temporal_angles/y_30_y1_300_y2_1000/";
    string save_file_folder = "/mnt/c/again_and_again/codes/git/plant_vertex_model/23_12_21_av_temporal/";
    int number_of_files_not_vtk =  9;
    std::vector<double> division_angles_all;
    std::vector<double> division_cell_y_all;
    std::vector<int> division_steps;
    std::vector<double> division_event_distance_to_tip_all;
    std::vector<double> max_y_per_step;
    std::vector<double> organ_length_per_step;
    //vector<double> minimum_curvature_vd;
    //vector<double> accumulated_negative_curvature_vd;
    //vector<double> tip_length_min_c_position_relative_vd;
    //vector<double> tip_length_min_c_position_absolute_vd;
    //std::string tip_length_min_c_position_relative_txt = save_file_folder+"tip_length_min_c_position_relative.txt";
    //std::string tip_length_min_c_position_absolute_txt = save_file_folder+"tip_length_min_c_position_absolute.txt";
    //std::string tip_length_min_c_position_relative_png = save_file_folder+"tip_length_min_c_position_relative.png";
    //std::string tip_length_min_c_position_absolute_png = save_file_folder+"tip_length_min_c_position_absolute.png";

    std::string divisionRecordFile = FileDirectory+"divisionRecord.txt";
        ifstream fin(divisionRecordFile,ios::in);
        if(!fin.is_open()){
            cout<<"Error: missing "<<divisionRecordFile <<" (divisionRecordFile!)"<<endl;        
            exit(1);
        }
        
        for(int i=0;i<10;i++){
            std::string tmp;
            fin>>tmp;
        }
        int division_event_number = 4941;
        int inner_cell_division_event=0;
        for(int i=0;i<division_event_number;i++){
            double double_tmp;
            fin>>double_tmp;
            int step_i;
            fin>>step_i;
            fin>>double_tmp;
            bool divided_cell_is_epidermal;
            fin>>divided_cell_is_epidermal;
            double division_angles_one;
            fin>>division_angles_one;
            fin>>double_tmp;
            double division_cell_y;
            fin>>division_cell_y;  
            fin>>double_tmp;
            long double long_double_tmp;
            fin>>long_double_tmp;
            fin>>long_double_tmp;
            if(divided_cell_is_epidermal==0){
                inner_cell_division_event++;
                division_steps.push_back(step_i/100000);
                division_angles_all.push_back(division_angles_one);
                division_cell_y_all.push_back(division_cell_y);
            }
        }
    
    int length_tmp = FileDirectory.length();
        char directory_tmp[length_tmp+1];
        strcpy(directory_tmp, FileDirectory.c_str());
        
        chdir(directory_tmp);
        int file_index = file_process::getFileNum(FileDirectory);
        cout<<"file_index: "<<file_index<<endl;
        if(file_index==0){
            cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
            exit(1);
        }
        cout<<"file_index "<<file_index<<endl;
        char str_cell[100], str_line[100];
        int simulation_steps = (file_index-number_of_files_not_vtk)/2;

    
        for(int i=0;i<simulation_steps;i++){
            Organ* p_g = new Organ;
            char str_cell[100], str_line[100];
            sprintf(str_cell,"2dv_cell%.5d00000.vtk",i);
            sprintf(str_line,"2dv_line%.5d00000.vtk",i);
            cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
            string CellVTK = FileDirectory+str_cell;
            string LineVTK = FileDirectory+str_line;
            readV::oneVTK(p_g, CellVTK, LineVTK,i);
            //int boundary_point_number = 100;
            //int boundary_point_neighbor=10;
            //int NumAveraging = 9;
            //geo::basic_VTK_analysis(p_g);
            //geo::calcGeometrics(p_g);
            organ_geo::organ_max_min_x_y(p_g);
            max_y_per_step.push_back(p_g->y_max_vertex);
            organ_length_per_step.push_back(p_g->y_max_vertex-p_g->y_min_vertex);
            std::cout<<"max_y_per_step "<<p_g->y_max_vertex<<"; organ length: "<<p_g->y_max_vertex-p_g->y_min_vertex<<std::endl;
            //vector<double> curvature_tmp = organ_geo::curvature_circle_fitting_kasa_along_polygon(p_g,boundary_point_number,boundary_point_neighbor,NumAveraging);            
            //int boundary_points_number = 100;
            //int fitting_distance = 10;
            //int NumAveraging = 9;
            //organ_geo::from_curvature_to_tip_position(p_g,boundary_points_number,fitting_distance,NumAveraging);
            //tip_length_min_c_position_relative_vd.push_back(p_g->tip_length_min_c_position_relative);
            //tip_length_min_c_position_absolute_vd.push_back(p_g->tip_length_min_c_position_absolute);

            //cout_fout_debug::fout_vector_double(tip_length_min_c_position_relative_vd,tip_length_min_c_position_relative_txt);
            //cout_fout_debug::fout_vector_double(tip_length_min_c_position_absolute_vd,tip_length_min_c_position_absolute_txt);

            //organ_geo::n_gons_analysis(p_g);
            /*
            ostringstream streamObj_1;
            streamObj_1 << fixed;
            streamObj_1 << setprecision(2);
            streamObj_1 << parameter_1[i];
            string parameter_1_str = streamObj_1.str();
            
            //cout<<"Pamater_1_str "<<parameter_1_str<<endl;
            //string save_curvature_txt = save_file_folder+string("_y_")+to_string(parameter_1[i])+string("_ab_")+to_string(parameter_2[j])+string("_")+to_string(k)+string(".txt");
            //string save_curvature_png = save_file_folder+string("_y_")+to_string(parameter_1[i])+string("_ab_")+to_string(parameter_2[j])+string("_")+to_string(k)+string(".png");
            string save_curvature_txt = save_file_folder+string("step_")+to_string(i)+string(".txt");
            string save_curvature_png = save_file_folder+string("step_")+to_string(i)+string(".png");
            cout<<"curvature output to "<<save_curvature_png<<endl;
            cout_fout_debug::fout_vector_double(curvature_tmp,save_curvature_txt);
            gnu_plot::curvature_plot(curvature_tmp,save_curvature_png);
            double minimum_curvature = geo_vv::vd_minimum(curvature_tmp);
            double accumulated_curvature = geo_vv::accumulated_negative(curvature_tmp);
            cout<<"minimum curvature "<<minimum_curvature<<"; accumulated negative curvature "<<accumulated_curvature<<endl;
            minimum_curvature_vd.push_back(minimum_curvature);
            accumulated_negative_curvature_vd.push_back(accumulated_curvature);
            */
            p_g->~Organ();
        }
        /*
        string minimum_curvature_txt = save_file_folder+"minimum_curvature.txt";
        string minimum_curvature_png = save_file_folder+"minimum_curvature.png";
        cout_fout_debug::fout_vector_double(minimum_curvature_vd,minimum_curvature_txt);
        string accumulated_negative_curvature_txt = save_file_folder + "accumulated_negative_curvature.txt";
        string accumulated_negative_curvature_png = save_file_folder + "accumulated_negative_curvature.png";
        cout_fout_debug::fout_vector_double(accumulated_negative_curvature_vd,accumulated_negative_curvature_txt);

        gnu_plot::curvature_plot(minimum_curvature_vd,minimum_curvature_png);
        gnu_plot::curvature_plot(accumulated_negative_curvature_vd,accumulated_negative_curvature_png);
        */

        

        for(int i=0;i<inner_cell_division_event;i++){
            for(int j=0;j<max_y_per_step.size();j++){
                if(division_steps[i]==j){
                    
                    double division_event_distance_to_tip = (max_y_per_step[j]-division_cell_y_all[i]);
                    division_event_distance_to_tip_all.push_back(division_event_distance_to_tip);
                    std::cout<<"For division event"<< i<<": Division raw y position "<<division_cell_y_all[i]<<" the maximum organ y position now "<<max_y_per_step[j]<<", the organ length now "<<organ_length_per_step[j];
                    std::cout<<"The normalized division position y "<<division_event_distance_to_tip<<std::endl;
                    break;
                }
            }
        }
        
        cout_fout_debug::fout_vector_double(division_angles_all,"division_angles_temporal.txt");
        cout_fout_debug::fout_vector_double(division_event_distance_to_tip_all,"division_event_distance_to_tip_temporal.txt");
}

else if(minor_mode=="single_time_lapse_basal_contour"){
    std::string input_file_path = "/mnt/e/c_file/2023_09_22_continue/y_30_y1_500_y2_500/";
    std::string input_geometrics_record = input_file_path + "geometrics_record.txt";

    int cell_number_min = 100;
    int cell_number_max = 16000;
    int cell_number_step = 300;
    double absolute_width_threshold = 28.0;
    vector<int> cell_numbers;
    for(int i=cell_number_min; i<=cell_number_max; i+=cell_number_step){
        cell_numbers.push_back(i);
    }

    std::vector<int> cell_number_per_step;
    std::vector<int> simulation_steps;
    std::ifstream file(input_geometrics_record);
    if (!file.is_open()) {
        std::cerr << "Unable to open file input_geometrics_record ("<<input_geometrics_record<<")." << std::endl;
        return 1;
    }
    std::string line;
    // Read file line by line
    while (getline(file, line)) {
        int x, y, z;

        std::istringstream iss(line);
        if (iss >> x >> y >> z) {  // Try to read the first two doubles from the line
            simulation_steps.push_back(x);
            cell_number_per_step.push_back(y+z);
            std::cout<<"step "<<x<<"; cell number: "<<y+z<<std::endl;
        } else {
            std::cerr << "Could not read two doubles from a line" << std::endl;
        }
    }
    file.close();

    std::vector<int> steps_for_cell_number;
    for(int i=0;i<cell_numbers.size();i++){
        int min_error=1000;
        int corresponding_step;
        for(int j=0;j<cell_number_per_step.size();j++){
            if(min_error>abs(cell_numbers[i]-cell_number_per_step[j]))
            {
             corresponding_step = simulation_steps[j];
             min_error=abs(cell_numbers[i]-cell_number_per_step[j]);
            }
        }
        steps_for_cell_number.push_back(corresponding_step);
        std::cout<<"For cell number "<<cell_numbers[i]<<" corresponding file step size is "<<steps_for_cell_number[i]<<std::endl;
    }

    std::vector<std::string> cell_vtk_files;
    std::vector<std::string> line_vtk_files;
    for(int i=0;i<cell_numbers.size();i++){
        char cell_vtk_file[100];
        sprintf(cell_vtk_file,"2dv_cell%.5d00000.vtk",steps_for_cell_number[i]/100000);
        cell_vtk_files.push_back(input_file_path+std::string(cell_vtk_file));
        char line_vtk_file[100];
        sprintf(line_vtk_file,"2dv_line%.5d00000.vtk",steps_for_cell_number[i]/100000);
        line_vtk_files.push_back(input_file_path+std::string(line_vtk_file));
        std::cout<<"For cell number "<<cell_numbers[i]<<" corresponding vtk files are "<<cell_vtk_files[i]<<","<<line_vtk_files[i]<<std::endl;
    }

    for(int step_i=0;step_i<cell_vtk_files.size();step_i++){

        char contour_txt_char[100], contour_png_char[100], contour_basal_txt_char[100], contour_basal_png_char[100], contour_absolute_width_txt_char[100],contour_absolute_width_png_char[100];
        sprintf(contour_txt_char,"contour_cell_number_%d.txt",cell_numbers[step_i]);
        sprintf(contour_png_char,"contour_cell_number_%d.png",cell_numbers[step_i]);
        sprintf(contour_basal_txt_char,"contour_basal_cell_number_%d.txt",cell_numbers[step_i]);
        sprintf(contour_basal_png_char,"contour_basal_cell_number_%d.png",cell_numbers[step_i]);
        sprintf(contour_absolute_width_txt_char,"contour_absolute_width_cell_number_%d.txt",cell_numbers[step_i]);
        sprintf(contour_absolute_width_png_char,"contour_absolute_width_cell_number_%d.png",cell_numbers[step_i]);
            Organ* p_g = new Organ;
            readV::oneVTK(p_g, cell_vtk_files[step_i], line_vtk_files[step_i]);

            geo::calcGeometrics(p_g);

            double sampling_distance = 0.001;
            double similarity_index;
            vector<Vertex> simulated_contour;
            for(int vi=0;vi<(int)p_g->p_v.size();vi++){
                if(p_g->p_v[vi]->IsSurface==1){
                    simulated_contour.push_back(*p_g->p_v[vi]);
                }
            }
            vector<Vertex> simulated_contour_normalized_y_only = geo_vv::normalization_y_only(simulated_contour);
            //vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
            vector<Vertex> simulated_outline_sampled_y_only = geo_vv::vector_vertex_sampling(simulated_contour_normalized_y_only,sampling_distance);
            //cout_fout_debug::cout_vector_vertex(simulated_outline_sampled_y_only);

            vector<Vertex> simulated_contour_normalized = geo_vv::normalization(simulated_contour);
            //vector<Vertex> simulated_contour_swapped = geo_vv::vv_x_swap(simulated_contour_normalized);
            vector<Vertex> simulated_outline_sampled = geo_vv::vector_vertex_sampling(simulated_contour_normalized,sampling_distance);

            std::vector<double> contour_widths;
            for(int width_i=0;width_i<1/sampling_distance;width_i++){
                double contour_width = abs(simulated_outline_sampled_y_only[width_i].loc.x) + abs(simulated_outline_sampled_y_only[2/sampling_distance-width_i].loc.x);
                contour_widths.push_back(contour_width);
            }
            cout_fout_debug::fout_vector_vertex(simulated_outline_sampled, std::string(contour_txt_char));
            gnu_plot::organ_contour_plot(std::string(contour_txt_char),std::string(contour_png_char),cell_numbers[step_i]);
            cout_fout_debug::fout_vector_double(contour_widths, std::string(contour_absolute_width_txt_char));
            gnu_plot::organ_width_plot(std::string(contour_absolute_width_txt_char),std::string(contour_absolute_width_png_char),cell_numbers[step_i]);
            vector<Vertex> basal_contour; //tip removed
            vector<int> tip_removed_index;
                //the very tip should be directly removed
                tip_removed_index.push_back(1/sampling_distance);


            for(int i=(int)1/sampling_distance-1; i>0; i-=sampling_distance){
                if(contour_widths[i]<absolute_width_threshold){
                    tip_removed_index.push_back(i);
                    tip_removed_index.push_back(2/sampling_distance-i);
                }
                else{
                    break;
                }
            }

            for(int i=0;i<simulated_outline_sampled.size();i++){
                bool is_true=1;
                for(int j=0;j<tip_removed_index.size();j++){
                    if(i==tip_removed_index[j]){
                        is_true=0;
                    }
                }

                if(is_true==1){
                    basal_contour.push_back(simulated_outline_sampled[i]);
                }
            }
            //cout_fout_debug::cout_vector_int(tip_removed_index);
            cout_fout_debug::fout_vector_vertex(basal_contour, std::string(contour_basal_txt_char));
            gnu_plot::organ_contour_plot(std::string(contour_basal_txt_char),std::string(contour_basal_png_char),cell_numbers[step_i]);
            p_g->~Organ();
    }
}

else if(minor_mode=="final"){

char current_directory[256];
getcwd(current_directory,256);

Organ *p_g = new Organ;

string FileDirectory = "/mnt/d/c_file/2023/2023_02_08_calculation_time_analysis/end_cell_10000/0/";
//string FileDirectory = "/mnt/d/c_file/2023/2023_02_23_biregion_angles_identity/y_0.60/y_0.60_ab_20.00/0/";

readV::final_VTK(p_g, FileDirectory);

geo::basic_VTK_analysis(p_g);

chdir(current_directory);

geo::calcGeometrics(p_g);


int boundary_points_number = 100;
int fitting_distance = 10;
int NumAveraging = 9;
organ_geo::from_curvature_to_tip_position(p_g,boundary_points_number,fitting_distance,NumAveraging);

}

else if(minor_mode=="batch_final_"){
cout<<"minor_mode: "<<minor_mode<<endl;

char current_directory[256];
getcwd(current_directory,256);

string FileDirectory = "/mnt/e/c_file/2023_06_20_biregion_angles_position_5000/";
string save_file_folder = "/mnt/c/again_and_again/codes/Git/analysis_result/2023_06_20_biregion_angles_position_5000/";
string analysis_filename="biregion_angles_position";
vector<string> plot_filenames;
//vector<double> parameter_1_plot = {0.80,0.60,0.40,0.20};
//vector<double> parameter_2_plot = {1.0,0.50,0,-0.50,-1.0};


//prepare to analyze the minimum curvature
int boundary_point_number=100;
int boundary_point_neighbor=10;
int NumAveraging = 9;

vector<double> circularity_batch;
vector<double> cell_number_ratio_batch;
vector<double> organ_area_batch;
vector<double> organ_circumference_batch;
vector<double> av_outermost_length_batch;
vector<double> av_epi_regularity_batch;
vector<double> av_in_regularity_batch;
vector<double> av_epi_in_area_ratio_batch;
vector<double> organ_length_batch;
vector<double> organ_width_batch;
vector<double> organ_leaf_index_batch;
vector<double> organ_cell_arrangement_x_batch;
vector<double> organ_cell_arrangement_y_batch;
vector<double> organ_cell_arrangement_ratio_batch;
vector<double> organ_potential_energy_batch;
vector<double> difference_index_batch;
vector<double> minimum_curvature_batch;
vector<double> accumulated_negative_curvature_batch;

//define batch parameters
int repeat_time =3;
string parameter_1_name="y";
string parameter_2_name="ab";


double parameter_1_start=0.00, parameter_1_end=1.00, parameter_1_step=0.10;
double parameter_2_start=0.05, parameter_2_end=1.0, parameter_2_step=0.05;


vector<double> parameter_1;
vector<double> parameter_2;
cout<<"parameter_1 "<<endl;
for(int i=parameter_1_start;i<((parameter_1_end-parameter_1_start)/parameter_1_step+1);i++){
    cout<<parameter_1_start+i*parameter_1_step<<endl;
    parameter_1.push_back(parameter_1_start+i*parameter_1_step);
}
//parameter_1.push_back(1.00);
cout<<"parameter_2 "<<endl;

for(int i=0;i<((parameter_2_end-parameter_2_start)/parameter_2_step+1);i++){
    //double parameter_2_tmp = pow(10,parameter_2_start+i*parameter_2_step);
    double parameter_2_tmp = parameter_2_start + i*parameter_2_step;
    cout<<parameter_2_tmp<<endl;
    //parameter_2.push_back(parameter_2_tmp);
}

parameter_2 = {0.00,0.10,0.20,0.30,0.40,0.50,0.75,1.00,1.50,2.00,3.00,4.00,5.00,7.50,10.00,20.00};
//parameter_2 = {0.00,0.10,0.20,0.30,0.40,0.50,0.75,1.00,2.00,3.00,4.00,5.00,7.50,10.00,20.00};
vector<double> parameter_1_plot = parameter_1;
vector<double> parameter_2_plot = parameter_2;

for(int i=0;i<parameter_1.size();i++){
    char parameter_1_file[100];
    sprintf(parameter_1_file,"y_%.2f/",parameter_1[i]);
    string parameter_1_file_s = parameter_1_file;
    for(int j=0;j<parameter_2.size();j++){
        char parameter_2_file[100];
        sprintf(parameter_2_file,"y_%.2f_ab_%.2f/",parameter_1[i],parameter_2[j]);
        string parameter_2_file_s = parameter_2_file;

        double circularity_av=0;
        double cell_number_ratio_av=0;
        double organ_area_av=0;
        double organ_circumference_av=0;
        double av_outermost_length_av=0;
        double av_epi_regularity_av=0;
        double av_in_regularity_av=0;
        double av_epi_in_area_ratio=0;
        double organ_length_av=0;
        double organ_width_av=0;
        double organ_leaf_index_av=0;
        double organ_cell_arrangement_x_av=0;
        double organ_cell_arrangement_y_av=0;
        double organ_cell_arrangement_ratio_av=0;
        double organ_potential_energy_av=0;
        double difference_index_av=0;
        double minimum_curvature_av=0;
        double accumulated_negative_curvature_av=0;
        

        sigma_O = 0.5;
        sigma_L = 2.0;
        kappa_S = 1.0;
        S_std=3.5;

        int zero_organ_width_count=0;
        for(int k=0;k<repeat_time;k++){
            cout<<"******************************************************"<<endl;
            char repeat_file[100];
            sprintf(repeat_file,"%d",k);
            string repeat_file_s = repeat_file;
            string directory_tmp2 = FileDirectory+parameter_1_file_s+parameter_2_file_s+repeat_file_s+"/";
            //PrintMemoryUsage();
            Organ *p_g = new Organ;
            cout<<"For "<<parameter_1_name<<"="<<parameter_1[i]<<", "<<parameter_2_name<<"="<<parameter_2[j]<<","<<"repeat time="<<k<<endl;
            //readV::final_VTK(p_g, directory_tmp2);
            
            string VTK_file_path = directory_tmp2;
            int length_tmp = VTK_file_path.length();
            char directory_tmp[length_tmp+1];
            strcpy(directory_tmp, VTK_file_path.c_str());
            
            chdir(directory_tmp);
            int file_index = file_process::getFileNum(VTK_file_path);
            cout<<"file_index: "<<file_index<<endl;
            if(file_index==0){
                cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
                exit(1);
            }
            //cout<<"file_index "<<file_index<<endl;
            char str_cell[100], str_line[100];
            sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-2);
            sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-2);
            cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
            string CellVTK = VTK_file_path+"/"+str_cell;
            string LineVTK = VTK_file_path+"/"+str_line;

            readV::oneVTK(p_g, CellVTK, LineVTK);



            geo::basic_VTK_analysis(p_g);

            geo::calcGeometrics(p_g);
            /*
            circularity_av+=p_g->circularity;
            cell_number_ratio_av +=(double)p_g->N_inner_cell/(double)p_g->N_epi_cell;
            organ_area_av+=p_g->area;
            
            organ_circumference_av+=p_g->perimeter;
            av_outermost_length_av+=p_g->perimeter/(double)p_g->N_epi_cell;
            av_epi_regularity_av+=p_g->regularity_epi_av;
            av_in_regularity_av+=p_g->regularity_in_av;


            organ_geo::organ_potential_energy(p_g);
            organ_potential_energy_av+=p_g->organ_potential_energy;
            cout<<"potential_energy: "<<p_g->organ_potential_energy<<endl;

            
            organ_geo::organ_leaf_index(p_g);
            organ_length_av+=p_g->organ_length;
            organ_geo::organ_cell_arrangement(p_g);
            if(p_g->organ_width==0){
                zero_organ_width_count++;
            }
            else{
                organ_width_av+=p_g->organ_width;
                organ_leaf_index_av+=p_g->leaf_index;
                organ_cell_arrangement_x_av+=p_g->cell_arrangement_x;
                organ_cell_arrangement_y_av+=p_g->cell_arrangement_y;
                organ_cell_arrangement_ratio_av+=p_g->cell_arrangement_ratio;
            }
            
            double av_epi_area, av_in_area=0;
            for(int ci=0;ci<(int)p_g->p_c.size();ci++){
                if(p_g->p_c[ci]->IsEpidermal==0){
                    av_in_area+=p_g->p_c[ci]->area;
                }
                else{
                    av_epi_area+=p_g->p_c[ci]->area;
                }
            }
            av_epi_area=av_epi_area/(double)p_g->N_epi_cell;
            av_in_area=av_in_area/(double)p_g->N_inner_cell;
            av_epi_in_area_ratio+=av_epi_area/av_in_area; 
            */
            vector<double> curvature_tmp = organ_geo::curvature_circle_fitting_kasa_along_polygon(p_g,boundary_point_number,boundary_point_neighbor,NumAveraging);

            
            //string save_curvature_txt = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_point_number)+string("_fitting_distance_")+to_string(boundary_point_neighbor)+string("_y_")+to_string(parameter_1[i])+string("_ab_")+to_string(parameter_2[j])+string("_")+to_string(k)+string(".txt");
            //string save_curvature_png = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_point_number)+string("_fitting_distance_")+to_string(boundary_point_neighbor)+string("_y_")+to_string(parameter_1[i])+string("_ab_")+to_string(parameter_2[j])+string("_")+to_string(k)+string(".png");
            ostringstream streamObj_1;
            streamObj_1 << fixed;
            streamObj_1 << setprecision(2);
            streamObj_1 << parameter_1[i];
            string parameter_1_str = streamObj_1.str();
            ostringstream streamObj_2;
            streamObj_2 << fixed;
            streamObj_2 << setprecision(2);
            streamObj_2 << parameter_2[j];
            string parameter_2_str = streamObj_2.str();
            
            cout<<"Pamater_1_str "<<parameter_1_str<<endl;
            //string save_curvature_txt = save_file_folder+string("_y_")+to_string(parameter_1_str[i])+string("_ab_")+to_string(parameter_2_str[j])+string("_")+to_string(k)+string(".txt");
            //string save_curvature_png = save_file_folder+string("_y_")+to_string(parameter_1_str[i])+string("_ab_")+to_string(parameter_2_str[j])+string("_")+to_string(k)+string(".png");
            string save_curvature_txt = save_file_folder+string("y_")+parameter_1_str+string("_ab_")+parameter_2_str+string("_")+to_string(k)+string(".txt");
            string save_curvature_png = save_file_folder+string("y_")+parameter_1_str+string("_ab_")+parameter_2_str+string("_")+to_string(k)+string(".png");

            cout<<"curvature output to "<<save_curvature_png<<endl;
            cout_fout_debug::fout_vector_double(curvature_tmp,save_curvature_txt);
            gnu_plot::curvature_plot(curvature_tmp,save_curvature_png);

            //if(k==0){
                string save_vtk_png = save_file_folder+string("y_")+parameter_1_str+string("_ab_")+parameter_2_str+string("_")+to_string(k)+string("_vtk.png");
                autoVTK::VTKLineCell(LineVTK, CellVTK,save_vtk_png);
            //}
            
            if(find(parameter_1_plot.begin(),parameter_1_plot.end(),parameter_1[i])!=parameter_1_plot.end()&&find(parameter_2_plot.begin(),parameter_2_plot.end(),parameter_2[j])!=parameter_2_plot.end()&&k==0){
                plot_filenames.push_back(save_curvature_png);    
            }

            
            double minimum_curvature = geo_vv::vd_minimum(curvature_tmp);
            minimum_curvature_av +=minimum_curvature;
            double accumulated_curvature = geo_vv::accumulated_negative(curvature_tmp);
            accumulated_negative_curvature_av += accumulated_curvature;
            cout<<"minimum curvature "<<minimum_curvature<<"; accumulated negative curvature "<<accumulated_curvature<<endl;




        }
        minimum_curvature_av/=repeat_time;
        minimum_curvature_batch.push_back(minimum_curvature_av);
        accumulated_negative_curvature_av/=repeat_time;
        accumulated_negative_curvature_batch.push_back(accumulated_negative_curvature_av);
        /*
        circularity_av/=repeat_time;
        circularity_batch.push_back(circularity_av);
        cell_number_ratio_av/=repeat_time;
        cell_number_ratio_batch.push_back(cell_number_ratio_av);
        organ_area_av/=repeat_time;
        organ_area_batch.push_back(organ_area_av);
        organ_circumference_av/=repeat_time;
        organ_circumference_batch.push_back(organ_circumference_av);
        av_outermost_length_av/=repeat_time;
        av_outermost_length_batch.push_back(av_outermost_length_av);
        av_epi_regularity_av/=repeat_time;
        av_epi_regularity_batch.push_back(av_epi_regularity_av);
        av_in_regularity_av/=repeat_time;
        av_in_regularity_batch.push_back(av_in_regularity_av);
        av_epi_in_area_ratio/=repeat_time;
        av_epi_in_area_ratio_batch.push_back(av_epi_in_area_ratio);
        organ_length_av/=repeat_time;
        organ_length_batch.push_back(organ_length_av);
        organ_width_av/=(repeat_time-zero_organ_width_count);
        organ_width_batch.push_back(organ_width_av);
        organ_leaf_index_av/=(repeat_time-zero_organ_width_count);
        organ_leaf_index_batch.push_back(organ_leaf_index_av);
        organ_cell_arrangement_x_av /=(repeat_time-zero_organ_width_count);
        organ_cell_arrangement_x_batch.push_back(organ_cell_arrangement_x_av);
        organ_cell_arrangement_y_av/=(repeat_time-zero_organ_width_count);
        organ_cell_arrangement_y_batch.push_back(organ_cell_arrangement_y_av);
        organ_cell_arrangement_ratio_av/=(repeat_time-zero_organ_width_count);
        organ_cell_arrangement_ratio_batch.push_back(organ_cell_arrangement_ratio_av);
        organ_potential_energy_av/=repeat_time;
        organ_potential_energy_batch.push_back(organ_potential_energy_av);
        */
    }
}

chdir(current_directory);


/*
string output_filename = analysis_filename+"_circularity.txt";
output::batch_final_analysis(parameter_1,parameter_2,circularity_batch,output_filename);
string output_filename2 = analysis_filename+"_cell_number_ratio.txt";
output::batch_final_analysis(parameter_1,parameter_2,cell_number_ratio_batch,output_filename2);
string output_filename3 = analysis_filename+"_organ_area.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_area_batch,output_filename3);
string output_filename4 = analysis_filename+"_organ_perimeter.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_circumference_batch,output_filename4);
string output_filename5 = analysis_filename+"_av_outermost_length.txt";
output::batch_final_analysis(parameter_1,parameter_2,av_outermost_length_batch,output_filename5);
string output_filename6 = analysis_filename+"_av_epi_regularity.txt";
output::batch_final_analysis(parameter_1,parameter_2,av_epi_regularity_batch,output_filename6);
string output_filename7 = analysis_filename+"_av_in_regularity.txt";
output::batch_final_analysis(parameter_1,parameter_2,av_in_regularity_batch,output_filename7);
string output_filename8 = analysis_filename+"_av_epi_in_area_ratio.txt";

output::batch_final_analysis(parameter_1,parameter_2,av_epi_in_area_ratio_batch,output_filename8);
string output_filename9 = analysis_filename+"_organ_length.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_length_batch,output_filename9);
string output_filename10 = analysis_filename+"_organ_width.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_width_batch,output_filename10);
string output_filename11 = analysis_filename+"_organ_leaf_index.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_leaf_index_batch,output_filename11);
string output_filename12 = analysis_filename+"_organ_cell_arrangement_x.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_cell_arrangement_x_batch,output_filename12);
string output_filename13 = analysis_filename+"_organ_cell_arrangement_y.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_cell_arrangement_y_batch,output_filename13);
string output_filename14 = analysis_filename+"_organ_cell_arrangement_ratio.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_cell_arrangement_ratio_batch,output_filename14);

string output_filename15 = analysis_filename+"_organ_potential_energy.txt";
output::batch_final_analysis(parameter_1,parameter_2,organ_potential_energy_batch,output_filename15);
*/

string output_filename16 = analysis_filename+"_organ_minimum_curvature.txt";
output::batch_final_analysis(parameter_1,parameter_2,minimum_curvature_batch,output_filename16);
string output_filename17 = analysis_filename+"_accumulated_negative_curvature.txt";
output::batch_final_analysis(parameter_1,parameter_2,accumulated_negative_curvature_batch,output_filename17);


//now the separated plotting is done, use openCV to automatically arrange them into a panel
    cout<<plot_filenames.size()<<" plots have done, now we are going to arrange them into a panel"<<endl;
    
    cout_fout_debug::cout_vector_string(plot_filenames);
        vector<cv::Mat> images;
        for (const auto& filename : plot_filenames) {
            images.push_back(cv::imread(filename));
        }
        
        // Create the panel image
        cv::Mat panel;
        
        // Arrange the images into a panel
        int images_per_row = parameter_1_plot.size();
        cout<<"images_per_row "<<images_per_row<<endl;
        int images_per_column = parameter_2_plot.size();
        cout<<"images_per_column "<<images_per_column<<endl;
        //int images_per_row = 7;
        //int images_per_column = 6;
        for (int i = 0; i < images_per_column; ++i) {
            cv::Mat row;
            for (int j = 0; j < images_per_row; ++j) {
                if (j == 0) {
                    row = images[i * images_per_row + j];
                } else {
                    cv::hconcat(row, images[i * images_per_row + j], row);
                }
            }
            if (i == 0) {
                panel = row;
            } else {
                cv::vconcat(panel, row, panel);
            }
        }
        cout<<"Panel is finished and we are going to save it !"<<endl;
        // Save the panel image
        string save_panel_name = save_file_folder+string("curvature_analysis_panel.png");
        cv::imwrite(save_panel_name, panel);

}

else if(minor_mode=="batch_final__"){
    std::string target_shape_contour = "/mnt/c/again_and_again/paper_writting/figure1/av_boundary_points_sampled.txt";
    std::vector<Vertex*> target_contour_sampled = boundary_geo::read_and_process_real_organ_contour_imagej(target_shape_contour);
    string FileDirectory = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/";
    string save_file_folder = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/analysis_result/contour/";
    std::vector<tuple<double,double,double>> heatmap_min_curvature_data;
    std::vector<tuple<double,double,double>> heatmap_similarity_data;
    std::vector<tuple<double,double,double>> heatmap_tip_length_relative_data;
    
    int repeat_time=3;
    std::string parameter_name_1 = "y";
    std::string parameter_name_2 = "db";

    std::vector<tuple<double,double,double>> tip_length_relative_all;
    std::vector<tuple<double,double,double>> min_curvature_all;
    std::vector<tuple<double,double,double>> similarity_all;
    

    vector<double> param_1 = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    //vector<double> param_2 = {10.0};
    //vector<double> param_1={0.5};
    vector<double> param_2 = {0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0};
    for(int i=0;i<param_1.size();i++){
        char parameter_1_file[100];
        sprintf(parameter_1_file,"y_%.2f/",param_1[i]);
        string parameter_1_file_s = parameter_1_file;
        for(int j=0;j<param_2.size();j++){
            char parameter_2_file[100];
            sprintf(parameter_2_file,"yb_%.2f_db_%.2f/",param_1[i],param_2[j]);
            string parameter_2_file_s = parameter_2_file;

            double tip_length_relative_av=0;
            double min_curvature_av=0;
            double similarity_av=0;
            for(int k=0;k<repeat_time;k++){
                cout<<"******************************************************"<<endl;
                char repeat_file[100];
                sprintf(repeat_file,"%d",k);
                string repeat_file_s = repeat_file;
                string directory_tmp2 = FileDirectory+parameter_1_file_s+parameter_2_file_s+repeat_file_s+"/";
                //PrintMemoryUsage();
                Organ *p_g = new Organ;
                cout<<"For "<<parameter_name_1<<"="<<param_1[i]<<", "<<parameter_name_2<<"="<<param_2[j]<<","<<"repeat time="<<k<<endl;

                string VTK_file_path = directory_tmp2;
                int length_tmp = VTK_file_path.length();
                char directory_tmp[length_tmp+1];
                strcpy(directory_tmp, VTK_file_path.c_str());
                
                chdir(directory_tmp);
                int file_index = file_process::getFileNum(VTK_file_path);
                cout<<"file_index: "<<file_index<<endl;
                if(file_index==0){
                    cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
                    exit(1);
                }
                //cout<<"file_index "<<file_index<<endl;
                char str_cell[100], str_line[100];
                sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-2);
                sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-2);
                cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
                string CellVTK = VTK_file_path+"/"+str_cell;
                string LineVTK = VTK_file_path+"/"+str_line;

                readV::oneVTK(p_g, CellVTK, LineVTK);
                std::cout<<"Cell number "<<p_g->p_c.size()<<std::endl;
                organ_geo::epidermal_identity(p_g);
                double similarity_index=boundary_geo::similarity_cal_during_simulation(p_g,target_contour_sampled);
                std::cout<<"Similarity_index is "<<similarity_index<<std::endl;
                similarity_all.push_back(std::make_tuple(param_1[i],param_2[j],similarity_index));
                similarity_av+=similarity_index/repeat_time;
                /*
                geo::calcGeometrics(p_g);
                int boundary_points_number = 100;
                int fitting_distance = 10;
                int NumAveraging = 9;
                organ_geo::from_curvature_to_tip_position(p_g,boundary_points_number,fitting_distance,NumAveraging);

                tip_length_relative_all.push_back(std::make_tuple(param_1[i],param_2[j],p_g->tip_length_min_c_position_relative));
                min_curvature_all.push_back(std::make_tuple(param_1[i],param_2[j],p_g->minimum_curvature));
                similarity_all.push_back(std::make_tuple(param_1[i],param_2[j],p_g->similarity_index));
                std::cout<<"minimum curvature "<<p_g->minimum_curvature<<std::endl;
                std::cout<<"tip length relative "<<p_g->tip_length_min_c_position_relative<<std::endl;
                std::cout<<"similarity index "<<p_g->similarity_index<<std::endl;
                tip_length_relative_av+=p_g->tip_length_min_c_position_relative/repeat_time;
                min_curvature_av+=p_g->minimum_curvature/repeat_time;
                similarity_av+=p_g->similarity_index/repeat_time;
                p_g->~Organ();

                char save_parameterc[100];
                sprintf(save_parameterc,"y_%.2f_db_%.2f_k_%d",param_1[i],param_2[j],k);
                std::string save_parameters = save_parameterc;
                autoVTK::VTKLineCell(LineVTK,CellVTK,save_file_folder+save_parameters+"_vtk.png");
                */
            }
            std::tuple<double,double,double> similarity_data = std::make_tuple(param_1[i],param_2[j],similarity_av);
            heatmap_similarity_data.push_back(similarity_data);
            /*
            std::tuple<double,double,double> min_curvature_data = std::make_tuple(param_1[i],param_2[j],min_curvature_av);
            std::tuple<double,double,double> tip_length_relative_data = std::make_tuple(param_1[i],param_2[j],tip_length_relative_av);
            
            heatmap_tip_length_relative_data.push_back(tip_length_relative_data);
            heatmap_min_curvature_data.push_back(min_curvature_data);
            
            */
        }
    }
    cout_fout_debug::fout_tuple_double(similarity_all,save_file_folder+"similarity_all.txt");
    cout_fout_debug::fout_tuple_double(heatmap_similarity_data,save_file_folder+"similarity_heatmap.txt");
    /*
    cout_fout_debug::fout_tuple_double(tip_length_relative_all,save_file_folder+"tip_length_relative_all.txt");
    cout_fout_debug::fout_tuple_double(min_curvature_all,save_file_folder+"min_curvature_all.txt");
    

    cout_fout_debug::fout_tuple_double(heatmap_min_curvature_data,save_file_folder+"min_curvature_heatmap.txt");
    cout_fout_debug::fout_tuple_double(heatmap_tip_length_relative_data,save_file_folder+"tip_length_relative_heatmap.txt");
    
    */
}

else if(minor_mode=="batch_final"){

char current_directory[256];
getcwd(current_directory,256);

string FileDirectory = "/mnt/e/c_file/2023_09_13_biregion_angles_gradual/";
string save_file_folder = "/mnt/e/c_file/2023_09_13_biregion_angles_gradual/";
string analysis_filename="y2_0.5";


//real_organ_contour_processed_for_similarity_index=boundary_geo::read_and_process_real_organ_contour_imagej(real_organ_contour_imagej_txt);

analysis_2D a_b;
//define batch parameters
std::cout<<"***********Define batch parameters******************"<<std::endl;

int repeat_time =1;
std::string parameter_name_1 = "y1";
std::string parameter_name_2 = "y2";

vector<double> param_1 = {100,200,300,500,1000,2000};
vector<double> param_2 = {100,200,300,500,1000,2000};

a_b.set_repeat_time(repeat_time);
a_b.add_param_name(parameter_name_1,parameter_name_2);
a_b.add_param_triang_right_up(param_1, param_2);
a_b.print();
a_b.add_all_file_directory_3(FileDirectory);
/*
string FileDirectory = "/mnt/e/c_file/2023_09_13_biregion_angles_gradual/";
string save_file_folder = "/mnt/c/again_and_again/codes/Git/analysis_result/2023_09_13_biregion_angles_gradual/";
string analysis_filename="biregion_angles_gradual";
vector<string> plot_filenames;

analysis_2D a_b;

//define batch parameters
std::cout<<"***********Define batch parameters******************"<<std::endl;

int repeat_time =3;
std::string parameter_name_1 = "y1";
std::string parameter_name_2 = "y2";
double param_1_init=0.4, param_1_term=1.0, param_1_incremental=0.1, param_2_init=0.4, param_2_term=1.0, param_2_incremental=0.1;

a_b.set_repeat_time(repeat_time);
a_b.add_param_name(parameter_name_1,parameter_name_2);
a_b.add_param_triang_right_up_incremental(param_1_init, param_1_term, param_1_incremental, param_2_init, param_2_term,param_2_incremental);
a_b.print();
a_b.add_all_file_directory(FileDirectory);
*/
for(int i=0;i<a_b.size();i++){
    cout<<"******************************************************"<<endl;
    std::cout<<"Analyzing parameters: "<<a_b.get_param_name().first<<"="<<a_b.get_parameter(i).first<<", "<<a_b.get_param_name().second<<"="<<a_b.get_parameter(i).second;
    std::cout<<" in "<<a_b.get_file_directory(i)<<std::endl;
    if(a_b.get_repeat_time()==1){
        std::string save_parameters = a_b.get_parameter_s(i);
        Organ *p_g = new Organ;
        std::pair<std::string,std::string> VTK_directory_tmp = a_b.VTK_directory(i);
        readV::oneVTK(p_g, VTK_directory_tmp);    
        autoVTK::VTKLineCell(VTK_directory_tmp,save_file_folder+save_parameters+"_vtk.png");
        geo::calcGeometrics_cout(p_g);        
        a_b.geometrics_recording(p_g,i,-1);
    }
    else{
        for(int k=0;k<a_b.get_repeat_time();k++){
            cout<<"******************************************************"<<endl;
            std::cout<<"For repeat time: "<<k<<std::endl;
            std::string save_parameters = a_b.get_parameter_s(i)+"_"+to_string(k);
            Organ *p_g = new Organ;
            std::pair<std::string,std::string> VTK_directory_tmp = a_b.VTK_directory(k,i);
            readV::oneVTK(p_g, VTK_directory_tmp);
            
            autoVTK::VTKLineCell(VTK_directory_tmp,save_file_folder+save_parameters+"_vtk.png");
            geo::calcGeometrics_cout(p_g);
            
            a_b.geometrics_recording(p_g,i,k);
        }
    }
    
    a_b.geometrics_output(save_file_folder);
}

}

else if(minor_mode=="batch_contours_smilarity"){
    cout<<"minor_mode: "<<minor_mode<<endl;

    string real_contour_file = "../analysis_example/mature_ts_contour_for_similarity_index.txt";

    string FileDirectory = "../../batch_contour/";
    string save_file_folder = "../../analysis_result/2023_06_29_batch_contour";
    string analysis_filename= "biregion_angles_position_";

    //string simulated_contour_file = "../../batch_contour/y_0.50_ab_5.00_0_outline.txt";

    //boundary_geo::similarity_Index_2(real_contour_file,simulated_contour_file,0.01);
    

    char current_directory[256];
    getcwd(current_directory,256);

    //define batch parameters
    int repeat_time =3;
    string parameter_1_name="y";
    string parameter_2_name="ab";


    vector<double> parameter_1 = {0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00};
    vector<double> parameter_2 = {0.00,0.10,0.20,0.30,0.40,0.50,0.75,1.00,1.50,2.00,3.00,4.00,5.00,7.50,10.00,20.00};
    //parameter_2 = {0.00,0.10,0.20,0.30,0.40,0.50,0.75,1.00,2.00,3.00,4.00,5.00,7.50,10.00,20.00};
    vector<double> parameter_1_plot = parameter_1;
    vector<double> parameter_2_plot = parameter_2;

    vector<double> similarity_index_batch;
    double best_similarity=1;
    double best_parameter_1;
    double best_parameter_2;
    double best_repeat_time;

    for(int i=0;i<parameter_1.size();i++){

        for(int j=0;j<parameter_2.size();j++){
            
            double similarity_index_av=0;
            for(int k=0;k<repeat_time;k++){
                cout<<"**************************************************************"<<endl;
                cout<<"For "<<parameter_1_name<<"="<<parameter_1[i]<<", "<<parameter_2_name<<"="<<parameter_2[j]<<", "<<"repeat time="<<k<<endl;;
                char simulated_contour_file[100];
                sprintf(simulated_contour_file,"y_%.2f_ab_%.2f_%d_outline.txt",parameter_1[i],parameter_2[j],k);
                string simulated_contour_file_s = simulated_contour_file;
                string directory_tmp2 = FileDirectory + simulated_contour_file_s;
                double similarity_index_tmp = boundary_geo::similarity_Index_2(real_contour_file,directory_tmp2,0.01);
                if(similarity_index_tmp<best_similarity){
                    best_similarity= similarity_index_tmp;
                    best_parameter_1 = parameter_1[i];
                    best_parameter_2 = parameter_2[j];
                    best_repeat_time = k;
                }
                similarity_index_av+=similarity_index_tmp;
            }
            similarity_index_av=similarity_index_av/(double)repeat_time;
            similarity_index_batch.push_back(similarity_index_av);
            cout<<"**************************************************************"<<endl;
            cout<<"For "<<parameter_1_name<<"="<<parameter_1[i]<<", "<<parameter_2_name<<"="<<parameter_2[j]<<", similarity index average="<<similarity_index_av<<endl;
        }
    }
    chdir(current_directory);
    string output_filename = analysis_filename+"similarity_index.txt";
    output::batch_final_analysis(parameter_1,parameter_2,similarity_index_batch,output_filename);
    cout<<"The best similarity is "<<best_similarity<<" and the parameter_1 is "<<best_parameter_1<<", parameter_2 is "<<best_parameter_2<<", repeat time is "<<best_repeat_time<<endl;
}

else if(minor_mode=="contour_extraction"){
    //std::string vtkFilePath = "/mnt/c/again_and_again/codes/Git/basic_simulations_results/circular_shape/";
    std::string vtkFilePath = "/mnt/e/c_file/2023_08_03_arrest_front_k/y_20_t_1200_t_-0.01/0/";
    std::string saveFilePath = "/mnt/c/again_and_again/codes/Git/analysis_result/2023_08_03_arrest_front_k_y_20_t_1200_t_-0.01/";
    int not_vtk_file_number = 2;
    int max_cell_number = 2500;


    int length_tmp = vtkFilePath.length();
        char directory_tmp[length_tmp+1];
        strcpy(directory_tmp, vtkFilePath.c_str());
        chdir(directory_tmp);
        int file_index = file_process::getFileNum(vtkFilePath);
        cout<<"file_index: "<<file_index<<endl;
        if(file_index==0){
            cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
            exit(1);
        }
        cout<<"file_index "<<file_index<<endl;
        char str_cell[100], str_line[100];
        int simulation_steps = (file_index-not_vtk_file_number)/2;
        cout<<"simulation_steps "<<simulation_steps<<endl;
        std::vector<int> cell_number;
        for(int i=0;i<simulation_steps;i++){
            Organ* p_g = new Organ;
            char str_cell[100], str_line[100];
            sprintf(str_cell,"2dv_cell%.5d00000.vtk",i);
            sprintf(str_line,"2dv_line%.5d00000.vtk",i);
            cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
            
            string CellVTK = vtkFilePath+str_cell;
            string LineVTK = vtkFilePath+str_line;
            readV::oneVTK(p_g, CellVTK, LineVTK,i);
            std::cout<<"Cell number is "<<p_g->p_c.size()<<std::endl;
            cell_number.push_back(p_g->p_c.size());
            organ_geo::epidermal_identity(p_g);
            char save_file_name[100];
            sprintf(save_file_name, "contour_%.5d00000.txt",i);
            std::string save_file_path = saveFilePath+save_file_name;
            output::simulated_contour(p_g,save_file_path);
            cout_fout_debug::fout_vector_int(cell_number,saveFilePath+"cell_number_per_step.txt");
            if(p_g->p_c.size()>max_cell_number){
                i=simulation_steps;
                break;
            }
            p_g->~Organ();
        }

        //cout_fout_debug::fout_vector_int(cell_number,"cell_number_per_step.txt");
}

else if(minor_mode=="batch_contour_extraction"){
string FileDirectory = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/";
    string save_file_folder = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/analysis_result/contour/";
    std::vector<tuple<double,double,double>> heatmap_min_curvature_data;
    std::vector<tuple<double,double,double>> heatmap_similarity_data;
    std::vector<tuple<double,double,double>> heatmap_tip_length_relative_data;
    
    int repeat_time=3;
    std::string parameter_name_1 = "y";
    std::string parameter_name_2 = "db";

    std::vector<tuple<double,double,double>> tip_length_relative_all;
    std::vector<tuple<double,double,double>> min_curvature_all;
    std::vector<tuple<double,double,double>> similarity_all;
    

    vector<double> param_1 = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    //vector<double> param_2 = {10.0};
    //vector<double> param_1={0.5};
    vector<double> param_2 = {0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0};
    for(int i=0;i<param_1.size();i++){
        char parameter_1_file[100];
        sprintf(parameter_1_file,"y_%.2f/",param_1[i]);
        string parameter_1_file_s = parameter_1_file;
        for(int j=0;j<param_2.size();j++){
            char parameter_2_file[100];
            sprintf(parameter_2_file,"yb_%.2f_db_%.2f/",param_1[i],param_2[j]);
            string parameter_2_file_s = parameter_2_file;

            double tip_length_relative_av=0;
            double min_curvature_av=0;
            double similarity_av=0;
            for(int k=0;k<repeat_time;k++){
                cout<<"******************************************************"<<endl;
                char repeat_file[100];
                sprintf(repeat_file,"%d",k);
                string repeat_file_s = repeat_file;
                string directory_tmp2 = FileDirectory+parameter_1_file_s+parameter_2_file_s+repeat_file_s+"/";
                //PrintMemoryUsage();
                Organ *p_g = new Organ;
                cout<<"For "<<parameter_name_1<<"="<<param_1[i]<<", "<<parameter_name_2<<"="<<param_2[j]<<","<<"repeat time="<<k<<endl;

                string VTK_file_path = directory_tmp2;
                int length_tmp = VTK_file_path.length();
                char directory_tmp[length_tmp+1];
                strcpy(directory_tmp, VTK_file_path.c_str());
                
                chdir(directory_tmp);
                int file_index = file_process::getFileNum(VTK_file_path);
                cout<<"file_index: "<<file_index<<endl;
                if(file_index==0){
                    cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
                    exit(1);
                }
                //cout<<"file_index "<<file_index<<endl;
                char str_cell[100], str_line[100];
                sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_index/2-2);
                sprintf(str_line,"2dv_line%.5d00000.vtk",file_index/2-2);
                cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
                string CellVTK = VTK_file_path+"/"+str_cell;
                string LineVTK = VTK_file_path+"/"+str_line;

                readV::oneVTK(p_g, CellVTK, LineVTK);
                std::cout<<"Cell number "<<p_g->p_c.size()<<std::endl;
                organ_geo::epidermal_identity(p_g);

                char str_save_contour[100], str_save_contour_png[100];
                sprintf(str_save_contour,"%s_%.2f_%s_%.2f_k_%d_contour.txt",parameter_name_1.c_str(),param_1[i],parameter_name_2.c_str(),param_2[j],k);
                sprintf(str_save_contour_png,"%s_%.2f_%s_%.2f_k_%d_contour.png",parameter_name_1.c_str(),param_1[i],parameter_name_2.c_str(),param_2[j],k);
                output::simulated_contour(p_g,save_file_folder+str_save_contour,save_file_folder+str_save_contour_png);
                p_g->~Organ();
            }
        }
    }
}

else if(minor_mode=="batch_simulated_simulated_similarity"){
    //from extracted contours and cell_number_per_step.txt
    string save_heatmap_name = "/mnt/e/2023_08_03_arrest_front_k/heatmap_k_-0.01_self_similarity.txt";

    std::string vtkPath1 = "/mnt/c/again_and_again/codes/Git/analysis_result/2023_08_03_arrest_front_k_y_20_t_1200_t_-0.01/";
    std::string vtkPath2 = "/mnt/c/again_and_again/codes/Git/analysis_result/2023_08_03_arrest_front_k_y_20_t_1200_t_-0.01/";
    std::string cell_number1_record = "/mnt/e/c_file/2023_08_03_arrest_front_k/y_20_t_1200_t_-0.01/0/GeometricsRecord.txt";
    std::string cell_number2_record = "/mnt/e/c_file/2023_08_03_arrest_front_k/y_20_t_1200_t_-0.01/0/GeometricsRecord.txt";

    int min_cell_number_1=0;
    int step_cell_number_1=100;
    int max_cell_number_1=2500;
    int min_cell_number_2=0;
    int step_cell_number_2=100;
    int max_cell_number_2=2500;

    std::vector<int> file1_cell_number;
    for(int i=min_cell_number_1;i<=max_cell_number_1;i+=step_cell_number_1){
        file1_cell_number.push_back(i);
    }

    std::vector<int> file2_cell_number;
    for(int i=min_cell_number_2;i<=max_cell_number_2;i+=step_cell_number_2){
        file2_cell_number.push_back(i);
    }

    std::vector<int> file1_step;
    std::vector<int> stepValues1;
    std::vector<int> cellnumberValues1;

    std::vector<int> file2_step;
    std::vector<int> stepValues2;
    std::vector<int> cellnumberValues2;
    std::string line;
    
    // Open file
    std::ifstream file(cell_number2_record);
    if (!file.is_open()) {
        std::cerr << "Unable to open file cell_number2_record ("<<cell_number2_record<<")." << std::endl;
        return 1;
    }
    // Read file line by line
    while (getline(file, line)) {
        int x, y, z;

        std::istringstream iss(line);
        if (iss >> x >> y >> z) {  // Try to read the first two doubles from the line
            stepValues2.push_back(x);
            cellnumberValues2.push_back(y+z);
            std::cout<<"step "<<x<<"; cell number: "<<y+z<<std::endl;
        } else {
            std::cerr << "Could not read two doubles from a line" << std::endl;
        }
    }
    file.close();

    for(int i=0;i<file2_cell_number.size();i++){
        int min_error=1000;
        int corresponding_step;
        for(int j=0;j<cellnumberValues2.size();j++){
            if(min_error>abs(file2_cell_number[i]-cellnumberValues2[j]))
            {
             corresponding_step = stepValues2[j];
             min_error=abs(file2_cell_number[i]-cellnumberValues2[j]);
            }
        }
        file2_step.push_back(corresponding_step);
        std::cout<<"For cell number "<<file2_cell_number[i]<<" corresponding file step size is "<<file2_step[i]<<std::endl;
    }

        // Open file
    std::ifstream file2(cell_number1_record);
    if (!file2.is_open()) {
        std::cerr << "Unable to open file cell_number1_record ("<<cell_number1_record<<")." << std::endl;
        return 1;
    }
    // Read file line by line
    while (getline(file2, line)) {
        int x, y, z;
        std::istringstream iss(line);
        if (iss >> x >> y >> z) {  // Try to read the first two doubles from the line
            stepValues1.push_back(x);
            cellnumberValues1.push_back(y+z);
            std::cout<<"step "<<x<<"; cell number: "<<y+z<<std::endl;
        } else {
            std::cerr << "Could not read two doubles from a line" << std::endl;
        }
    }

    file2.close();

    for(int i=0;i<file1_cell_number.size();i++){
        int min_error=1000;
        int corresponding_step;
        for(int j=0;j<cellnumberValues1.size();j++){
            if(min_error>abs(file1_cell_number[i]-cellnumberValues1[j]))
            {
             corresponding_step = stepValues1[j];
             min_error=abs(file1_cell_number[i]-cellnumberValues1[j]);
            }
        }
        file1_step.push_back(corresponding_step);
        std::cout<<"For cell number "<<file1_cell_number[i]<<" corresponding file step size is "<<file1_step[i]<<std::endl;
    }
    //std::vector<int> file2_cell_number = {200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4200,4400,4600,,4800,5000,5200,5400,5600,5800,6000,6200,6400,6600,6800,7000,7200,7400,7600,7800,8000,8200,8400,8600,8800,9000};

    //std::vector<int> file2_step = {17,30,41,54,64,74,81,88,93,99,108,117,126,134,141,160,178,197,216,236,262};
    //std::vector<int> file2_cell_number = {100,200,300,400,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,4500,5000};
    
    //std::vector<int> file2_step = {19,31,39,46,54,67,77,92,104,124,142,165,181,196,211,226,240};
    //std::vector<int> file2_cell_number = {100,200,300,400,500,750,1000,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000};
    std::vector<double> difference_index;

    std::vector<tuple<int,int,double>> heatmap_data;
    for(int i=0;i<file1_step.size();i++){
        char file1_char[100];
        sprintf(file1_char,"contour_%.5d00000.txt",file1_step[i]);
        for(int j=0;j<file2_step.size();j++){
            char file2_char[100];
            sprintf(file2_char,"contour_%.5d00000.txt",file2_step[j]);    
            std::string file1 = vtkPath1+file1_char;
            std::string file2 = vtkPath2+file2_char;
            double similarity_tmp = boundary_geo::similarity_for_contour_contour(file1,file2,0.01);
            std::cout<<"For file "<<file1<<" and "<<file2<<" the difference index is "<< similarity_tmp<<std::endl;
            difference_index.push_back(similarity_tmp);
            heatmap_data.push_back(make_tuple(file1_cell_number[i],file2_cell_number[j],similarity_tmp));
        }
    }
    //cout_fout_debug::fout_heatmap_data(file1_cell_number,file2_cell_number,difference_index,"heatmap.txt");
    
    cout_fout_debug::fout_tuple_double(heatmap_data,save_heatmap_name);
    
    /*
    FILE *pipe = popen("gnuplot -persist", "w"); // open a pipe to Gnuplot

    if (pipe == nullptr) {
        std::cerr << "Could not open pipe" << std::endl;
        return 1;
    }
    fprintf(pipe, "set terminal png\n"); // set the terminal to PNG
    fprintf(pipe, "set output 'heatmap.png'\n"); // set the output file name
    fprintf(pipe, "set palette rgb 35,13,10\n");
    //fprintf(pipe, "set view map\n"); // set the Gnuplot to a heatmap view
    fprintf(pipe, "plot '-' with image notitle\n"); // plot data

    for (size_t i = 0; i < file1_cell_number.size(); ++i) {
        fprintf(pipe, "%d %d %lf\n", file1_cell_number[i], file2_cell_number[i], difference_index[i]); // send the data to Gnuplot
    }

    fprintf(pipe, "e\n"); // end of data
    fflush(pipe); // flush the pipe

    std::cout << "Heatmap displayed, press enter to exit." << std::endl;
    std::cin.get(); // wait for the user to press enter

    pclose(pipe); // close the pipe
    
*/

    //vector<double> minimum_curvature_vd;
    //vector<double> accumulated_negative_curvature_vd;
    
}

else if(minor_mode=="time_laspe_basal_contour_similarity"){
    //double similarity= boundary_geo::similarity_for_contour_contour("/mnt/c/again_and_again/codes/git/plant_vertex_model/y_30_basal_part_similarity/contour_basal_cell_number_2800.txt","/mnt/c/again_and_again/codes/git/plant_vertex_model/y_40_basal_part_similarity/contour_basal_cell_number_5500.txt",0.001);
    
    
    string save_heatmap_txt = "/mnt/c/again_and_again/codes/git/plant_vertex_model/y_40_basal_part_similarity/y_20_30_heatmap.txt";

    std::string vtkPath1 = "/mnt/c/again_and_again/codes/git/plant_vertex_model/y_20_basal_part_similarity/";
    std::string vtkPath2 = "/mnt/c/again_and_again/codes/git/plant_vertex_model/y_30_basal_part_similarity/";
    double sampling_distance = 0.001;

    int min_cell_number_1=1000;
    int step_cell_number_1=300;
    int max_cell_number_1=10000;
    int min_cell_number_2=1000;
    int step_cell_number_2=300;
    int max_cell_number_2=16000;

    std::vector<int> cell_numbers1;
    std::vector<int> cell_numbers2;
    std::vector<std::string> vtkFiles1;
    std::vector<std::string> vtkFiles2;
    for(int i=min_cell_number_1;i<=max_cell_number_1;i+=step_cell_number_1){
        char vtkFile1[100];
        sprintf(vtkFile1,"contour_basal_cell_number_%d.txt",i);
        vtkFiles1.push_back(vtkPath1+std::string(vtkFile1));
        cell_numbers1.push_back(i);
    }
    for(int i=min_cell_number_2;i<=max_cell_number_2;i+=step_cell_number_2){
        char vtkFile2[100];
        sprintf(vtkFile2,"contour_basal_cell_number_%d.txt",i);
        vtkFiles2.push_back(vtkPath2+std::string(vtkFile2));
        cell_numbers2.push_back(i);
    }
    cout_fout_debug::cout_vector_string(vtkFiles1);
    cout_fout_debug::cout_vector_string(vtkFiles2);
    cout_fout_debug::cout_vector_int(cell_numbers1);
    cout_fout_debug::cout_vector_int(cell_numbers2);
    std::vector<tuple<int,int,double>> heatmap_data;
    std::vector<double> similarities;
    for(int i=0;i<vtkFiles1.size();i++){
        for(int j=0;j<vtkFiles2.size(); j++){
            double similarity= boundary_geo::similarity_for_contour_contour(vtkFiles1[i],vtkFiles2[j],sampling_distance);
            similarities.push_back(similarity);
            heatmap_data.push_back(make_tuple(cell_numbers1[i],cell_numbers2[j],similarity));
            std::cout<<"For "<<vtkFiles1[i]<<" and "<<vtkFiles2[j]<<" the difference index is "<<similarity<<std::endl;
        }
    }

        //cout_fout_debug::fout_heatmap_data(cell_numbers1,cell_numbers2,similarities,save_heatmap_txt);
        cout_fout_debug::fout_tuple_double(heatmap_data,save_heatmap_txt);
        //gnu_plot::heatmap(save_heatmap_txt,save_heatmap_png,cell_numbers[step_i]);
        
}

else if(minor_mode=="repeat_time_lapse_contour_similarity_heatmap_video"){
    std::string filePath = "/mnt/c/again_and_again/codes/git/2023_10_22_temporal_TS_repeat_analysis/";
    int repeat_time = 10;


    for(int i=0;i<repeat_time;i++){
        std::string filePathi = filePath + std::to_string(i) + std::string("/");
        
    }
}

else if(minor_mode=="geometrics_txt"){
    std::string geometrics_txt_path = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/y10/0/geometrics_record.txt";
    std::string save_file_path = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/y10_";
    Geometrics_analysis_class ga = readV::read_geo_output_txt(geometrics_txt_path);
    
    std::string figure_x_axis = "cell_number";
    std::string figure_y_axis = "length_width_ratio";

    std::string save_file_name = save_file_path + "_x_" + figure_x_axis + "_y_" + figure_y_axis + ".txt";
    geo_analysis::geo_plot(ga, figure_x_axis, figure_y_axis, save_file_name);
    
}

else if(minor_mode=="geometrics_txt_batch"){
    std::string FileDirectory = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/";
    std::string save_file_path = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/";
    
    std::string parameter_name = "y";
    int repeat_time=3;
    std::vector<int> param = {10,11,12,13,14,15,16,17,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,34,36,38,40};
    std::vector<int> param_;
    std::vector<double> transition_time_all;
    std::vector<double> transition_time_av;
    for(size_t i=0;i<param.size();i++){
        char parameter_file[100];
        sprintf(parameter_file,"y%d/",param[i]);
        std::string parameter_file_s = parameter_file; 
        double av_transition_time=0;
        for(int k=0;k<repeat_time;k++){
            char repeat_file[100];
            sprintf(repeat_file,"%d",k);
            string repeat_file_s=repeat_file;
            std::string directory_tmp = FileDirectory+parameter_file_s+repeat_file_s+"/";
            std::string geometrics_txt_path = directory_tmp+"geometrics_record.txt";
            std::cout<<"For file "<<parameter_name<<"="<<param[i]<<",repeat_time="<<k<<std::endl;
            Geometrics_analysis_class ga=readV::read_geo_output_txt(geometrics_txt_path);
            std::cout<<"This txt should have "<<ga.value.size()-1<<"lines "<<std::endl;
            std::vector<double> length_width_ratio_tmp = ga.getColumnData("length_width_ratio");
            int transition_time;
            for(int step_i=0;step_i<ga.value.size()-1;step_i++){
                if(areAllValuesGreaterThan(length_width_ratio_tmp,step_i,1.1)){
                    transition_time=step_i;
                    break;
                }
            }

            int transition_time_cell_number = ga.getColumnData("inner_cell_number")[transition_time]+ga.getColumnData("epidermal_cell_number")[transition_time];
            transition_time_all.push_back(transition_time_cell_number);
            param_.push_back(param[i]);
            av_transition_time+=transition_time_cell_number/(double)repeat_time;
            std::cout<<"The transition time by cell number is "<<transition_time_cell_number<<std::endl;
        }
        
        transition_time_av.push_back(av_transition_time);
        std::cout<<"The av transition time is "<<av_transition_time<<std::endl;
    }

    cout_fout_debug::fout_vector_int_double(param_,transition_time_all,save_file_path+"transition_time_all.txt");
    cout_fout_debug::fout_vector_int_double(param,transition_time_av,save_file_path+"transition_time_av.txt");
    /*
    std::string parameter_name_1 = "y";
    std::string parameter_name_2 = "db";
    int repeat_time = 3;
    vector<double> param_1 = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    vector<double> param_2 = {0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.5,10.0};

    //std::vector<double> heatmap_x;
    //std::vector<double> heatmap_y;
    std::vector<double> similarity_index_vec;
    
    for(int i=0;i<param_1.size();i++){
        char parameter_1_file[100];
        sprintf(parameter_1_file,"y_%.2f/",param_1[i]);
        string parameter_1_file_s = parameter_1_file;
        for(int j=0;j<param_2.size();j++){
            char parameter_2_file[100];
            sprintf(parameter_2_file,"yb_%.2f_db_%.2f/",param_1[i],param_2[j]);
            string parameter_2_file_s = parameter_2_file;

            std::vector<double> similarity_index;
            for(int k=0;k<repeat_time;k++){
                cout<<"******************************************************"<<endl;
                char repeat_file[100];
                sprintf(repeat_file,"%d",k);
                string repeat_file_s = repeat_file;
                string directory_tmp2 = FileDirectory+parameter_1_file_s+parameter_2_file_s+repeat_file_s+"/";
                cout<<"For "<<parameter_name_1<<"="<<param_1[i]<<", "<<parameter_name_2<<"="<<param_2[j]<<","<<"repeat time="<<k<<endl;

                std::string geometrics_txt_path = directory_tmp2 + "geometrics_record.txt";
                Geometrics_analysis_class ga = readV::read_geo_output_txt(geometrics_txt_path);
                std::cout<<"This txt should have "<<ga.value.size()-1<<"lines "<<std::endl;
                std::cout<<"Similarity index is "<<ga.getColumnData("similarity_index")[ga.value.size()-1]<<std::endl;
                similarity_index.push_back(ga.getColumnData("similarity_index")[ga.value.size()-1]);
            }
            similarity_index_vec.push_back(wVec::v_av(similarity_index));
            std::cout<<"av similarity index is "<<wVec::v_av(similarity_index)<<std::endl;
        }
    }

    std::string save_file_name = save_file_path + "similarity_index_heatmap.txt";
    cout_fout_debug::fout_heatmap_data(param_2,param_1,similarity_index_vec,save_file_name);
    */
}

else if(minor_mode=="repeat_time_lapse_contour"){
    std::string input_file_path = "/mnt/e/c_file/2023_10_22_temporal_TS_repeat_analysis/";
    std::string save_file_path = "/mnt/e/c_file/2023_10_22_temporal_TS_repeat_analysis/";

    int repeat_time = 10;
    int cell_number_min = 0;
    int cell_number_step = 300;
    int cell_number_max = 15000;
    std::vector<int> cell_numbers;
    for(int i=cell_number_min;i<=cell_number_max;i+=cell_number_step){
        cell_numbers.push_back(i);
    }
    cout_fout_debug::cout_vector_int(cell_numbers);
    std::vector<std::vector<std::string>> contour_file_all;
    
    
    for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
        char path_i[100];
        sprintf(path_i,"%d/",repeat_i);
        std::string input_file_path_i = input_file_path + std::string(path_i);
        std::string input_cell_number_file = input_file_path_i + "/cell_number_per_step.txt";
        std::vector<std::string> contour_files;
        std::vector<int> contour_steps;

        //read from cell_number_per_step.txt
        std::vector<double> input_cell_numbers;
        std::vector<double> input_simulation_steps;
        std::ifstream file(input_cell_number_file);
        if (!file.is_open()) {
            std::cerr << "Unable to open file cell_number_per_step.txt ("<<input_cell_number_file<<")." << std::endl;
            return 1;
        }
        std::string line;
        file>>line;
        file>>line;
        double value1,value2;
        // Read file line by line
        while (file >> value1 >> value2) {
            input_simulation_steps.push_back(value1);
            input_cell_numbers.push_back(value2);
        }
        file.close();
        for(int i=0;i<cell_numbers.size();i++){
            int min_error=1000;
            int corresponding_step;
            for(int j=0;j<input_cell_numbers.size();j++){
                if(min_error>abs(cell_numbers[i]-input_cell_numbers[j]))
                {
                corresponding_step = input_simulation_steps[j];
                min_error=abs(cell_numbers[i]-input_cell_numbers[j]);
                }
            }
            contour_steps.push_back(corresponding_step);
            std::cout<<"For cell number "<<cell_numbers[i]<<" corresponding file step size is "<<contour_steps[i]<<std::endl;
        }

        for(int i=0;i<cell_numbers.size();i++){
            char contour_file[100];
            sprintf(contour_file,"contour_%.5d00000.txt",contour_steps[i]);
            contour_files.push_back(input_file_path_i+std::string(contour_file));
            std::cout<<"For cell number "<<cell_numbers[i]<<" corresponding contour file is "<<contour_files[i]<<std::endl;
        }

        contour_file_all.push_back(contour_files);
    }
    
    for(int step_i=0; step_i<cell_numbers.size();step_i++){
        std::vector<std::string> contour_files_per_step;
        for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
            contour_files_per_step.push_back(contour_file_all[repeat_i][step_i]);
        }
        char save_png_char[100];
        sprintf(save_png_char,"contours_cell_number_%d.png",cell_numbers[step_i]);
        std::string save_png = save_file_path + std::string(save_png_char);
        gnu_plot::multiple_time_lapse_contour_txt_plot(contour_files_per_step,save_png,cell_numbers[step_i]);
    }

    for(int step_i=0; step_i<cell_numbers.size();step_i++){
        std::vector<std::string> contour_files_per_step;
        for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
            contour_files_per_step.push_back(contour_file_all[repeat_i][step_i]);
        }
        std::vector<int> repeats;
        for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
            repeats.push_back(repeat_i);
        }
        std::vector<double> similarities;
        char save_heatmap_txt_char[100];
        sprintf(save_heatmap_txt_char,"heatmap_cell_number_%d.txt",cell_numbers[step_i]);
        std::string save_heatmap_txt = save_file_path + std::string(save_heatmap_txt_char);
        char save_heatmap_png_char[100];
        sprintf(save_heatmap_png_char,"heatmap_cell_number_%d.png",cell_numbers[step_i]);
        std::string save_heatmap_png = save_file_path + std::string(save_heatmap_png_char);

        for(int repeat_i=0;repeat_i<repeat_time;repeat_i++){
            for(int repeat_j=0;repeat_j<repeat_time;repeat_j++){
                double similarity_tmp   = boundary_geo::similarity_for_contour_contour(contour_files_per_step[repeat_i],contour_files_per_step[repeat_j],0.01);
                //std::cout<<"Difference index between "<<contour_files_per_step[repeat_i]<<","<<contour_files_per_step[repeat_j]<<" is "<<similarity_tmp<<std::endl;
                similarities.push_back(similarity_tmp);
                
            }
        }

        //cout_fout_debug::cout_vector_double(similarities);
        cout_fout_debug::fout_heatmap_data(repeats,repeats,similarities,save_heatmap_txt);
        gnu_plot::heatmap(save_heatmap_txt,save_heatmap_png,cell_numbers[step_i]);
    }

}

else if(minor_mode=="width_change"){
    std::string input_file_path = "/mnt/c/again_and_again/codes/git/plant_vertex_model/basal_part_similarity/";
    std::string input_file = input_file_path + "contour_width_10000_cell.txt";

    std::string output_file = input_file_path + "contour_width_change_10000_cell.txt";
    std::vector<double> widths;
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Unable to open file ("<<input_file<<")." << std::endl;
        return 1;
    }
    std::string line;
    file>>line;
    file>>line;
    double value1,value2;
    // Read file line by line
    while (file >> value1 >> value2) {
        widths.push_back(value2);
    }
    file.close();

    std::vector<double> width_changes;
    int averaged_steps = 10;

    for(int i=0;i<widths.size()-averaged_steps;i++){
        double width_change = (widths[i+averaged_steps] - widths[i])/averaged_steps;
        width_changes.push_back(width_change);
    }

    cout_fout_debug::fout_vector_double(width_changes,output_file);
}

else{
    cout<<"Fatal error: no minor mode selected ! (major mode is Analysis)"<<endl;
    exit(-1);
}

}

else if(major_mode=="experiment"){

if(minor_mode=="contour_simple"){
    string xy_outline_file = "/mnt/c/again_and_again/paper_writting/Fig9S_concave_apex/O.triangularis_raw_contour.txt";
    string save_file_folder = "/mnt/c/again_and_again/paper_writting/Fig9S_concave_apex/O.";
    int boundary_points_number = 100;
    int fitting_distance = 10;
    int averaing_points = 9;

    vector<Vertex*> xy_outline_vertex_imageJ = readV::xy_txt_to_vertexXY(xy_outline_file);
    cout_fout_debug::cout_vector_vertex(xy_outline_vertex_imageJ);
    vector<Vertex*> preprocessed_outline = geo_vv::after_ImageJ_process(xy_outline_vertex_imageJ);

        vector<Vertex*> contour_equal_distance_points = geo_vv::organ_boundary_points_along_polygon(preprocessed_outline,boundary_points_number);
        string save_boundary_equal_distance_txt = save_file_folder+string("contour_equal_distance_")+to_string(boundary_points_number)+string(".txt");
        string save_boundary_equal_distance_png = save_file_folder+string("contour_equal_distance_")+to_string(boundary_points_number)+string(".png");
        cout_fout_debug::fout_vector_vertex(contour_equal_distance_points,save_boundary_equal_distance_txt);
        gnu_plot::organ_contour_plot(contour_equal_distance_points,1,save_boundary_equal_distance_png);

        vector<double> curvature_ij = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(contour_equal_distance_points,fitting_distance,3);
        string save_curvature_txt = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_points_number)+string("_fitting_distance_")+to_string(fitting_distance)+string(".txt");
        string save_curvature_png = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_points_number)+string("_fitting_distance_")+to_string(fitting_distance)+string(".png");
        cout_fout_debug::fout_vector_double(curvature_ij,save_curvature_txt);
        gnu_plot::curvature_plot(curvature_ij,save_curvature_png);
}

else if(minor_mode=="averaged_contour"){
    std::string input_filefolder = "/mnt/c/again_and_again/paper_writting/figure1/contour/";
    std::vector<std::string> contour_files = wSystem::files_inside_folder(input_filefolder);
        //wVec::v_print(contour_files);
    std::vector<std::vector<Vertex>> contours;
    std::vector<std::vector<Vertex>> normalized_contours;
    std::vector<std::vector<Vertex>> sampled_contours;
    for(auto contour_file : contour_files){
        std::cout<<"For contour "<<contour_file<<std::endl;
    //0. read contour data from txt file and because it comes from imageJ, we have to do a vertical reflection
        std::vector<Vertex> contour = vVertex::read_from_txt(contour_file,1,2);
        contour = vVertex::vertical_reflection(contour);
        contours.push_back(contour);
            //vVertex::print(contour);

    //1. get normalized contour 
        std::vector<Vertex> normalized_contour = vVertex::normalization(contour);
            //vVertex::print(normalized_contour);
        normalized_contours.push_back(normalized_contour);
    //2. get the sampled contour
        std::vector<Vertex> sampled_contour = vVertex::sample(normalized_contour);    
        sampled_contours.push_back(sampled_contour);
            //gnuplot::vv(sampled_contour);
            //vVertex::print(sampled_contour);
    }


    //3. get the averaged contour from the sampled contour
    std::vector<Vertex> averaged_contour = vVertex::vec_averaged(sampled_contours);
            //vVertex::print(averaged_contour);
        //gnuplot::vv(averaged_contour);
    
    //4. before calculate curvature, we have extract the boundary points (with equal distance)


    //5. calculate the curvature and output it
    std::vector<Vertex> averaged_boundary_points =geo_vv::organ_boundary_points_along_polygon(averaged_contour,100);
    vector<double> curvature_averaged = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(averaged_boundary_points,10,3);
    //gnuplot::curvature_plot(curvature_averaged);
    std::cout<<"minimum curvature "<<geo_vv::vd_minimum(curvature_averaged)<<std::endl;
    vVertex::save_txt("/mnt/c/again_and_again/paper_writting/figure1/av_boundary_points.txt",averaged_boundary_points);
    cout_fout_debug::fout_vector_double(curvature_averaged,"/mnt/c/again_and_again/paper_writting/figure1/curvature_averaged.txt");
}

else if(minor_mode=="background_remove_single"){
        std::string input_image = "/mnt/c/again_and_again/processed_experiment_data/2024_07_24_oxalis_debilis_triangular/240724-OT-1-cut-1.tif";
        std::string output_image = "/mnt/c/again_and_again/processed_experiment_data/2024_07_24_oxalis_debilis_triangular/240724-OT-1-cut-1_processed.tif";
        cv::Mat image = cv::imread(input_image);

        if(image.empty()){
            std::cout << "Could not open or find the image" << std::endl;
            return -1;
        }

        // Define the range of grey color you want to replace
        cv::Scalar greyMin(0, 0, 0); // Adjust these values according to your image's grey
        cv::Scalar greyMax(55, 55, 55); // Adjust these values according to your image's grey

        // Create a mask that captures the areas of the image that are grey
        cv::Mat mask;
        cv::inRange(image, greyMin, greyMax, mask);

        // Replace the grey color in the original image with black
        image.setTo(cv::Scalar(0, 0, 0), mask);

        cv::Mat dst;
        Image_2dv::removeSmallerObjects(image, dst);
        // Save the resultant image
        cv::imwrite(output_image, dst);
}

else if(minor_mode=="background_remove_batch"){

    std::string input_file_folder = "/mnt/c/again_and_again/paper_writting/figure1/scaled_image/";
    std::string output_file_folder = "/mnt/c/again_and_again/paper_writting/figure1/background_removed/";
    for(const auto& entry:std::filesystem::directory_iterator(input_file_folder)){
        auto path=entry.path();
        std::string input_image = path.string();
        std::string output_image = output_file_folder+path.filename().string();
        cv::Mat image = cv::imread(path);

        if(image.empty()){
            std::cout << "Could not open or find the image" << std::endl;
            return -1;
        }

        // Define the range of grey color you want to replace
        cv::Scalar greyMin(0, 0, 0); // Adjust these values according to your image's grey
        cv::Scalar greyMax(90, 90, 90); // Adjust these values according to your image's grey

        // Create a mask that captures the areas of the image that are grey
        cv::Mat mask;
        cv::inRange(image, greyMin, greyMax, mask);

        // Replace the grey color in the original image with black
        image.setTo(cv::Scalar(0, 0, 0), mask);

        // Save the resultant image
        cv::imwrite(output_image, image);
    }
}

else if(minor_mode=="contour"){

//batch curvature analysis for a sinlge contour: testing parameters
    string xy_outline_file = "/mnt/c/again_and_again/paper_writting/figure1/1_0_contour.txt";
    string save_file_folder = "/mnt/c/again_and_again/paper_writting/figure1/";
    vector<string> plot_filenames;

    
    vector<double> boundary_points_number = {500,400,300,200,150,100,50};
    vector<double> fitting_distance = {1,2,5,10,20,30};

    //1. read and preprocess of xy txt file => preprocessed_outline

    vector<Vertex*> xy_outline_vertex_imageJ = readV::xy_txt_to_vertexXY(xy_outline_file);
    vector<Vertex*> preprocessed_outline = geo_vv::after_ImageJ_process(xy_outline_vertex_imageJ);

    for(int i=0;i<boundary_points_number.size();i++){
        int boundary_points_number_i = boundary_points_number[i];

    //2. calculate the contour made by sampled equal-distance points, and output the vertex position into txt files, then plot it 
        
        vector<Vertex*> contour_equal_distance_points = geo_vv::organ_boundary_points_along_polygon(preprocessed_outline,boundary_points_number_i);
        string save_boundary_equal_distance_txt = save_file_folder+string("contour_equal_distance_")+to_string(boundary_points_number_i)+string(".txt");
        string save_boundary_equal_distance_png = save_file_folder+string("contour_equal_distance_")+to_string(boundary_points_number_i)+string(".png");
        cout_fout_debug::fout_vector_vertex(contour_equal_distance_points,save_boundary_equal_distance_txt);
        gnu_plot::organ_contour_plot(contour_equal_distance_points,1,save_boundary_equal_distance_png);
        plot_filenames.push_back(save_boundary_equal_distance_png);
    //3. calculate the curvature of the contour points, and output the curvature into txt files, then plot it 
    
        for(int j=0;j<fitting_distance.size();j++){
            int fitting_distance_j = fitting_distance[j];
            vector<double> curvature_ij = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(contour_equal_distance_points,fitting_distance_j,1);
            string save_curvature_txt = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_points_number_i)+string("_fitting_distance_")+to_string(fitting_distance_j)+string(".txt");
            string save_curvature_png = save_file_folder+string("curvature_boundary_points_")+to_string(boundary_points_number_i)+string("_fitting_distance_")+to_string(fitting_distance_j)+string(".png");
            cout_fout_debug::fout_vector_double(curvature_ij,save_curvature_txt);
            gnu_plot::curvature_plot(curvature_ij,save_curvature_png);
            plot_filenames.push_back(save_curvature_png);
        }
    }

    //4. now the separated plotting is done, use openCV to automatically arrange them into a panel
    cout<<plot_filenames.size()<<" plots have done, now we are going to arrange them into a panel"<<endl;
    cout_fout_debug::cout_vector_string(plot_filenames);
        vector<cv::Mat> images;
        for (const auto& filename : plot_filenames) {
            images.push_back(cv::imread(filename));
        }
        
        // Create the panel image
        cv::Mat panel;
        
        // Arrange the images into a panel
        int images_per_row = fitting_distance.size()+1;
        int images_per_column = boundary_points_number.size();
        //int images_per_row = 7;
        //int images_per_column = 6;
        for (int i = 0; i < images_per_column; ++i) {
            cv::Mat row;
            for (int j = 0; j < images_per_row; ++j) {
                if (j == 0) {
                    row = images[i * images_per_row + j];
                } else {
                    cv::hconcat(row, images[i * images_per_row + j], row);
                }
            }
            if (i == 0) {
                panel = row;
            } else {
                cv::vconcat(panel, row, panel);
            }
        }

        // Save the panel image
        string save_panel_name = save_file_folder+string("curvature_analysis_panel.png");
        cv::imwrite(save_panel_name, panel);
    
    


//single curvature analysis 
/*
    int boundary_points_number = 200;
    int fitting_distance=10;
    int NumAveraging = 3;
    string curvature_output_file = save_file_folder+ "curvature.png";
    //using outline x_y_coordinates txt extracted from imageJ
    //1. read and preprocess of xy txt file

    vector<Vertex*> xy_outline_vertex_imageJ = readV::xy_txt_to_vertexXY(xy_outline_file);
    //cout_fout_debug::cout_vector_vertex(xy_outline_vertex_imageJ);
    vector<Vertex*> xy_outline_vertex_processed = geo_vv::after_ImageJ_process(xy_outline_vertex_imageJ);
    //cout_fout_debug::fout_vector_vertex(xy_outline_vertex_processed,"../analysis_result/circle_JW_normalized.txt");

    //calculate the organ area and organ perimeter from the xy_outline
    double vv_area = geo_vv::area_vv_boundary(xy_outline_vertex_processed);
    double vv_perimeter = geo_vv::perimeter_vv_boundary(xy_outline_vertex_processed);
    cout<<"area: "<<vv_area<<"; perimeter: "<<vv_perimeter<<endl;
    double vv_geometric_complexity = 0.25*vv_perimeter/sqrt(vv_area);
    cout<<"geometric entropy "<<vv_geometric_complexity<<endl;
    vector<Vertex*> xy_equal_distance_points = geo_vv::organ_boundary_points_along_polygon(xy_outline_vertex_processed,boundary_points_number);
    cout_fout_debug::cout_vector_vertex(xy_equal_distance_points);
    cout_fout_debug::fout_vector_vertex(xy_equal_distance_points,"../analysis_result/arabidopsis_serration_100_points.txt");
    vector<double> curvature_vv = geo_vv::curvature_circle_fitting_kasa_three_boundary_points(xy_equal_distance_points,fitting_distance,NumAveraging);
    //cout_fout_debug::fout_vector_double(curvature_vv,"../analysis_result/circle_JW_curvature.txt");
    gnu_plot::curvature_plot(curvature_vv,curvature_output_file);
    double minimum_curvature = geo_vv::vd_minimum(curvature_vv);
    double accumulated_curvature = geo_vv::accumulated_negative(curvature_vv);
    cout<<"minimum_curvature: "<<minimum_curvature<<"; accumulated_curvature: "<<accumulated_curvature<<endl;
*/

}

else if(minor_mode=="EdU"){
    string EdU_pair_file = "/mnt/d/UT_confocal_data/2023_06_05_TS_mature/TS_primordia_1_Cycle/EdU_pairs.txt";
    string image_name = "/mnt/d/UT_confocal_data/2023_06_05_TS_mature/TS_primordia_1_Cycle/angles_distribution_10.jpg";
    string angles_file = "/mnt/d/UT_confocal_data/2023_06_05_TS_mature/TS_primordia_1_Cycle/angles_10.txt";
    //1. read the csv file
    vector<Vertex*> EdU_points = readV::xy_txt_to_vertexXY(EdU_pair_file);

    //2. calculate leaf axis
    double leaf_axis_radians = atan2((EdU_points[1]->loc-EdU_points[0]->loc).y,(EdU_points[1]->loc-EdU_points[0]->loc).x);
    //change radians to degree
    double leaf_axis_degrees = leaf_axis_radians*(180.0/M_PI);
    cout<<"Leaf axis angles (degrees): "<<leaf_axis_degrees<<endl;

    //3. calculate the raw division direction
    vector<double> EdU_pairs(EdU_points.size()/2-1,0.0);
    for(int i=0;i<EdU_pairs.size();i++){
        EdU_pairs[i] = atan2((EdU_points[2*i+3]->loc-EdU_points[2*i+2]->loc).y,(EdU_points[2*i+3]->loc-EdU_points[2*i+2]->loc).x);
        EdU_pairs[i] = EdU_pairs[i]*(180.0/M_PI);
        if(EdU_pairs[i]<0){
            EdU_pairs[i]+=180;
        }
    }
    cout_fout_debug::cout_vector_double(EdU_pairs);
    //4. adjust the division direction by leaf axis
    vector<double> EdU_angles(EdU_pairs.size(),0.0);
    for(int i=0;i<EdU_angles.size();i++){
        EdU_angles[i] = EdU_pairs[i]-leaf_axis_degrees+90;
        EdU_angles[i] = abs(EdU_angles[i]-90); 
    }
    cout_fout_debug::cout_vector_double(EdU_angles);
    //5. plot a histogram of cell division angles distribution
    int binAngle = 10;
    double minAngle = 0.0;
    double maxAngle = 90.0;
    double binCount = (maxAngle-minAngle)/binAngle;
    vector<int> binAngles(binCount,0);
    for(double angle : EdU_angles){
        int binIndex = static_cast<int>((angle-minAngle)/binAngle);
        ++binAngles[binIndex];
    }
    cout_fout_debug::cout_vector_int(binAngles);
        FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe != NULL) {
        // Set the terminal to png and output to a file
        fprintf(pipe, "set terminal png\n");
        fprintf(pipe, "set output '%s'\n",image_name.c_str());
        fprintf(pipe, "set style fill solid border lc rgb 'black'\n");

        // Plot with points
        fprintf(pipe, "plot '-' with boxes lw 2 lc rgb 'light-cyan' notitle \n");
        // Output data points
        for(int i = 0; i < binAngles.size(); i++) {
            int binMiddle = minAngle + i*binAngle + binAngle/2;
            fprintf(pipe, "%d %d\n", binMiddle, binAngles[i]);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
    }
    //6. save the analyzed results
    cout_fout_debug::fout_vector_double(EdU_angles,angles_file);
}

else if(minor_mode=="EdU_axis"){
    std::string filefolder = "/mnt/e/UT_confocal_data/2023_08_16_EdU_TS/";
    string left_pairs_file = filefolder+"apical_pairs.txt";
    string right_pairs_file =  filefolder+"basal_pairs.txt";
    string left_axis_file =  filefolder+"axis.txt";
    string right_axis_file =  filefolder+"axis.txt";

    string image_name =  filefolder+"apical.png";
    //string angles_file = "/mnt/c/again_and_again/codes/Git/plant_vertex_model/analysis_example/2022_11_29_3_EdU/EdU_angles_all.txt";

    string left_norm_anlges_file =  filefolder+"apical_normalized_angles.txt";
    string right_norm_anlges_file =  filefolder+"basal_normalized_angles.txt";

    //read csv files
    vector<Vertex*> left_pairs = readV::xy_txt_to_vertexXY(left_pairs_file);
    //vector<Vertex*> right_pairs = readV::xy_txt_to_vertexXY(right_pairs_file);
    vector<Vertex*> left_axis = readV::xy_txt_to_vertexXY(left_axis_file);
    //vector<Vertex*> right_axis = readV::xy_txt_to_vertexXY(right_axis_file);

    vector<double> left_norm_angles = line_geo::angles_normalization_changing_axis(left_pairs, left_axis);
    //vector<double> right_norm_angles = line_geo::angles_normalization_changing_axis(right_pairs, right_axis);

    cout_fout_debug::fout_vector_double(left_norm_angles,left_norm_anlges_file);
    //cout_fout_debug::fout_vector_double(right_norm_angles,right_norm_anlges_file);

    /*
    int binAngle = 15;
    double minAngle = 0.0;
    double maxAngle = 180.0;
    double binCount = (maxAngle-minAngle)/binAngle;
    vector<int> binAngles(binCount,0);
    for(double angle : left_norm_angles){
        int binIndex = static_cast<int>((angle-minAngle)/binAngle);
        ++binAngles[binIndex];
    }
    for(double angle : right_norm_angles){
        int binIndex = static_cast<int>((angle-minAngle)/binAngle);
        ++binAngles[binIndex];
    }
    cout_fout_debug::cout_vector_int(binAngles);
        FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe != NULL) {
        // Set the terminal to png and output to a file
        fprintf(pipe, "set terminal png\n");
        fprintf(pipe, "set output '%s'\n",image_name.c_str());
        fprintf(pipe, "set style fill solid border lc rgb 'black'\n");
        fprintf(pipe, "set yrange [0:5]\n");
        fprintf(pipe, "set xtics '%d' \n",binAngle);

        // Plot with points
        fprintf(pipe, "plot '-' with boxes lw 2 lc rgb 'light-cyan' notitle \n");
        // Output data points
        for(int i = 0; i < binAngles.size(); i++) {
            int binMiddle = minAngle + i*binAngle + binAngle/2;
            fprintf(pipe, "%d %d\n", binMiddle, binAngles[i]);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
        
}*/
}

else if(minor_mode=="real_leaf_contour_averaged_contour_difference_index_comparison"){
    std::string input_filefolder = "/mnt/c/again_and_again/paper_writting/figure1/contour/";
    std::vector<std::string> contour_files = wSystem::files_inside_folder(input_filefolder);
        //wVec::v_print(contour_files);
    std::vector<std::vector<Vertex>> contours;
    std::vector<std::vector<Vertex>> normalized_contours;
    std::vector<std::vector<Vertex>> sampled_contours;
    for(auto contour_file : contour_files){
        std::cout<<"For contour "<<contour_file<<std::endl;
    //0. read contour data from txt file and because it comes from imageJ, we have to do a vertical reflection
        std::vector<Vertex> contour = vVertex::read_from_txt(contour_file,1,2);
        contour = vVertex::vertical_reflection(contour);
        contours.push_back(contour);
            //vVertex::print(contour);

    //1. get normalized contour 
        std::vector<Vertex> normalized_contour = vVertex::normalization(contour);
            //vVertex::print(normalized_contour);
        normalized_contours.push_back(normalized_contour);
    //2. get the sampled contour
        std::vector<Vertex> sampled_contour = vVertex::sample(normalized_contour);    
        sampled_contours.push_back(sampled_contour);
            //gnuplot::vv(sampled_contour);
            //vVertex::print(sampled_contour);
    //3. calculate the difference heatmap for real leaves
    }
    std::vector<double> similarity_index_vec;

    //averaged contour 
    std::vector<Vertex> averaged_contour = vVertex::vec_averaged(sampled_contours);
    for(int i=0;i<contour_files.size();i++){
        //heatmap_x.push_back(i);
        double similarity_index_tmp = boundary_geo::similarity_Index_1(sampled_contours[i],averaged_contour,0.01);
        cout<<"For leaf "<<contour_files[i]<<" the difference index with averaged contour is "<<similarity_index_tmp<<std::endl;
        similarity_index_vec.push_back(similarity_index_tmp);
    }
    cout_fout_debug::fout_vector_double(similarity_index_vec,"/mnt/c/again_and_again/paper_writting/figure1/difference_index_with_av_contours.txt");
}

else if(minor_mode=="real_leaf_contour_difference_index_comparison"){
    std::string input_filefolder = "/mnt/c/again_and_again/paper_writting/figure1/contour/";
    std::vector<std::string> contour_files = wSystem::files_inside_folder(input_filefolder);
        //wVec::v_print(contour_files);
    std::vector<std::vector<Vertex>> contours;
    std::vector<std::vector<Vertex>> normalized_contours;
    std::vector<std::vector<Vertex>> sampled_contours;
    for(auto contour_file : contour_files){
        std::cout<<"For contour "<<contour_file<<std::endl;
    //0. read contour data from txt file and because it comes from imageJ, we have to do a vertical reflection
        std::vector<Vertex> contour = vVertex::read_from_txt(contour_file,1,2);
        contour = vVertex::vertical_reflection(contour);
        contours.push_back(contour);
            //vVertex::print(contour);

    //1. get normalized contour 
        std::vector<Vertex> normalized_contour = vVertex::normalization(contour);
            //vVertex::print(normalized_contour);
        normalized_contours.push_back(normalized_contour);
    //2. get the sampled contour
        std::vector<Vertex> sampled_contour = vVertex::sample(normalized_contour);    
        sampled_contours.push_back(sampled_contour);
            //gnuplot::vv(sampled_contour);
            //vVertex::print(sampled_contour);
    //3. calculate the difference heatmap for real leaves
    }
    std::vector<double> similarity_index_vec;
    std::vector<int> heatmap_x;
    for(int i=0;i<contour_files.size();i++){
        heatmap_x.push_back(i);
        for(int j=0;j<contour_files.size();j++){
            double similarity_index_tmp = boundary_geo::similarity_Index_1(sampled_contours[i],sampled_contours[j],0.01);
            cout<<"For leaf "<<contour_files[i]<<" and "<<contour_files[j]<<" the difference index is "<<similarity_index_tmp<<std::endl;
            similarity_index_vec.push_back(similarity_index_tmp);
        }
    }

    cout_fout_debug::fout_heatmap_data(heatmap_x,heatmap_x,similarity_index_vec,"/mnt/c/again_and_again/paper_writting/figure1/similarity_heatmap.txt");
    
}

else if(minor_mode=="cell_shape___"){
    std::string file_folder = "/mnt/c/again_and_again/paper_writting/Fig3_biased_expansion/trial/";
    std::vector<std::string> files = {"3-1"};
    for(auto file:files){
        //std::string file = "1-1-epi";
        std::string file_extension = ".jpg";
        std::string cell_shape_file = file_folder+file+file_extension;
        std::cout<<"For file "<<file<<std::endl;
        cv::Mat img = cv::imread(cell_shape_file,cv::IMREAD_GRAYSCALE);
        if(img.empty()){
            std::cerr<<"Could not read the image: "<<cell_shape_file<<std::endl;
            exit(1);
        }
        //std::cout<<"pixel value for origin 0,0 "<<static_cast<int>(img.at<uchar>(0,0));
        //std::cout<<"pixel value for origin 1,0 "<<static_cast<int>(img.at<uchar>(0,1));
        //std::cout<<"pixel value for origin 400,300 "<<static_cast<int>(img.at<uchar>(400,300));
        int top, bottom, left, right;
        top = bottom = 500; // This will add 50 pixels to the top and bottom
        left = right = 500;  // This will add 50 pixels to the left and right

        // Create a new image with the border
        cv::Mat imageWithBorder;
        cv::copyMakeBorder(img, imageWithBorder, top, bottom, left, right, cv::BORDER_REPLICATE);

        cv::Mat edges;
        cv::Canny(imageWithBorder, edges, 20, 220);
        cv::namedWindow("Edges", cv::WINDOW_AUTOSIZE);
        cv::imshow("Edges", edges);
        cv::waitKey(0);

        

        // Find contours
        std::vector<std::vector<cv::Point>> contours;
        cv::findContours(edges, contours, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE);
        std::cout<<"Number of cells "<<contours.size()<<std::endl;
        double area = cv::contourArea(contours[0]);
            //std::cout<<"For cell "<<0<< "its area is "<<area<<"; ";
        //std::cout<<"Cell area: "<<cell_area<<std::endl;
        std::vector<double> cell_area;
        std::vector<double> length_width_ratio;
        std::vector<double> length_axis;
        // Draw contours on a blank image
        //cv::Mat contourImage = cv::Mat::zeros(img.size(), CV_8UC3);

        cv::Mat contourImage = cv::Mat::zeros(img.size(), CV_8UC3);
                    //cv::Scalar color = cv::Scalar(255, 255, 255); // White color
                    //cv::drawContours(contourImage, contours[0], static_cast<int>(0), color, 1, cv::LINE_8);
                    //cv::namedWindow("Contours", cv::WINDOW_AUTOSIZE);
                    //cv::imshow("Contours", contourImage);
                    //cv::waitKey(0);
        
        for (size_t i = 0; i < contours.size(); i++) {
            cv::Scalar color = cv::Scalar(255, 255, 255); // White color
            cv::drawContours(contourImage, contours, static_cast<int>(i), color, 1, cv::LINE_8);
            double area = cv::contourArea(contours[i]);
            std::cout<<"For cell "<<i<< "its area is "<<area<<"; ";
            if(area>200&&area<15000){
                
                cell_area.push_back(area);

                //cv::Rect boundingBox = cv::boundingRect(contours[i]);

                // Calculate the ratio of the width to the height of the bounding box
                //double ratio = static_cast<double>(boundingBox.height) / boundingBox.width;
                //std::cout<<"length/width ratio "<<ratio<<"; ";
                //length_width_ratio.push_back(ratio);
                cv::RotatedRect minRect = cv::minAreaRect(contours[i]);

                // The angle returned by minAreaRect is the angle between the longest side
                // of the rectangle and the horizontal axis. It is in the range [-90, 0].
                // When the angle is 0, the rectangle is horizontal.
                // When the angle is -90, the rectangle is vertical.
                float angle = minRect.angle;
                double ratio = minRect.size.width / minRect.size.height;

                if (minRect.size.width < minRect.size.height) {
                    
                    angle = 90 + angle; // Adjust the angle if the height is greater than the width
                    ratio = minRect.size.height / minRect.size.width;
                }
                angle = angle +90;
                if(angle>180){
                    angle = angle-180;
                }
                length_axis.push_back(angle);
                
                std::cout<<"length/width ratio "<<ratio<<"; ";
                length_width_ratio.push_back(ratio);
                std::cout<<" angles toward x-axis"<<angle<<std::endl;
            }
        }

        // Display the contour image
        //cv::namedWindow("Contours", cv::WINDOW_AUTOSIZE);
        //cv::imshow("Contours", contourImage);
        //cv::waitKey(0);
        std::string cell_analysis_txt = file_folder+file+".txt";
        std::string cell_area_png = file_folder+file+"_cell_area.png";
        std::string cell_biased_expansion_png = file_folder+file+"_biased_expansion.png";
        cout_fout_debug::fout_vector_double(cell_area,length_width_ratio,length_axis,cell_analysis_txt,"cell_area","length_width_ratio","length_axis");
        gnu_plot::histogram(cell_area,cell_area_png);
        gnuplot::vv(length_axis,length_width_ratio,cell_biased_expansion_png);
        
    }

}

else if(minor_mode=="cell_shape"){
    std::string file_folder = "/mnt/c/again_and_again/paper_writting/Fig3_biased_expansion/raw data/mature_epi_api/";
    //std::vector<std::string> files = {"1-epi_skeleton","2-epi_skeleton","3-epi_skeleton","4-epi_skeleton","5-epi_skeleton","6-epi_skeleton","1-pali_skeleton","2-pali_skeleton","3-pali_skeleton","4-pali_skeleton","5-pali_skeleton","6-pali_skeleton"};
    std::vector<std::string> files = {"2-1"};
    std::string cell_analysis_all_txt = file_folder+"cp_pali.txt";
    std::vector<double> cell_area_all;
    std::vector<double> length_width_ratio_all;
    std::vector<double> length_axis_all;
    for(auto file:files){
        //std::string file = "1-1-epi";
        std::string file_extension = ".jpg";
        std::string cell_shape_file = file_folder+file+file_extension;
        std::cout<<"For file "<<file<<std::endl;
        Image2dv* image = new Image2dv;
        Image_2dv::read_image(cell_shape_file,image);
        //Image_2dv::draw_image_from_points(image);
        Image_2dv::find_polygon_from_image(image);
        //Image_2dv::find_contours_from_polygon(image);
        //Image_2dv::statistics_of_polygons(image);
    }
}



else if(minor_mode=="cell_shape_"){
    std::string file_folder = "/mnt/c/again_and_again/paper_writting/Fig3_biased_expansion/raw data/mature_epi_api/";
    //std::vector<std::string> files = {"1-epi_skeleton","2-epi_skeleton","3-epi_skeleton","4-epi_skeleton","5-epi_skeleton","6-epi_skeleton","1-pali_skeleton","2-pali_skeleton","3-pali_skeleton","4-pali_skeleton","5-pali_skeleton","6-pali_skeleton"};
    std::vector<std::string> files = {"1-1"};
    //std::vector<std::string> files = {"1-3","2-3","2-4","2-5","2-6","3-3","3-4","3-5","3-6"};
    std::string cell_analysis_all_txt = file_folder+"cp_pali.txt";
    std::vector<double> cell_area_all;
    std::vector<double> length_width_ratio_all;
    std::vector<double> length_axis_all;
    for(auto file:files){
        //std::string file = "1-1-epi";
        std::string file_extension = ".jpg";
        std::string cell_shape_file = file_folder+file+file_extension;
        std::cout<<"For file "<<file<<std::endl;
        cv::Mat img = cv::imread(cell_shape_file,cv::IMREAD_GRAYSCALE);
        if(img.empty()){
            std::cerr<<"Could not read the image: "<<cell_shape_file<<std::endl;
            exit(1);
        }
        cv::Mat edges;
        cv::Canny(img, edges, 50, 150);
        cv::namedWindow("Edges", cv::WINDOW_AUTOSIZE);
        //cv::imshow("Edges", edges);
        //cv::waitKey(0);

        std::vector<double> cell_area;
        std::vector<double> length_width_ratio;
        std::vector<double> length_axis;

        // Find contours
        std::vector<std::vector<cv::Point>> contours;
        cv::findContours(edges, contours, cv::RETR_LIST, cv::CHAIN_APPROX_SIMPLE);
        std::cout<<"Number of cells "<<contours.size()<<std::endl;
        // Draw contours on a blank image
        
        for (size_t i = 0; i < contours.size(); i++) {

            double area = cv::contourArea(contours[i]);
            //if(area>300&&area<20000){
                if(contours[i].front() == contours[i].back()){
                    std::cout<<"Fatal error: the contour is not closed!"<<std::endl;
                }
                std::cout<<"For cell "<<i<< "its area is "<<area<<"; ";
                cell_area.push_back(area);
                    
                //cv::RotatedRect minRect = cv::minAreaRect(contours[i]);
                //std::vector<cv::Point> approx;
                //double epsilon = 0.01 * cv::arcLength(contours[i], true);
                //cv::approxPolyDP(contours[i], approx, epsilon, true);
                cv::RotatedRect minRect = cv::minAreaRect(contours[i]);
                float angle = minRect.angle;
                double ratio = minRect.size.width / minRect.size.height;
                if (minRect.size.width < minRect.size.height) {
                    
                    angle = 90 + angle; // Adjust the angle if the height is greater than the width
                    ratio = minRect.size.height / minRect.size.width;
                }
                else{
                    //angle+=90;
                }
                length_axis.push_back(angle);
                
                std::cout<<"length/width ratio "<<ratio<<"; ";
                length_width_ratio.push_back(ratio);
                std::cout<<" angles toward x-axis"<<angle<<std::endl;

                cell_area_all.push_back(area);
                length_axis_all.push_back(angle);
                length_width_ratio_all.push_back(ratio);
                    cv::Mat contourImage = cv::Mat::zeros(img.size(), CV_8UC3);
                    cv::Scalar color = cv::Scalar(255, 255, 255); // White color
                    cv::drawContours(contourImage, contours, static_cast<int>(i), color, 1, cv::LINE_8);
                    cv::namedWindow("Contours", cv::WINDOW_AUTOSIZE);
                    cv::imshow("Contours", contourImage);
                    cv::waitKey(0);
            //}
        }

        // Display the contour image
        
        std::string cell_analysis_txt = file_folder+file+".txt";
        std::string cell_area_png = file_folder+file+"_cell_area.png";
        std::string cell_biased_expansion_png = file_folder+file+"_biased_expansion.png";
        cout_fout_debug::fout_vector_double(cell_area,length_width_ratio,length_axis,cell_analysis_txt,"cell_area","length_width_ratio","length_axis");
        gnu_plot::histogram(cell_area,cell_area_png);
        gnuplot::vv(length_axis,length_width_ratio,cell_biased_expansion_png);
    }
    
    cout_fout_debug::fout_vector_double(cell_area_all,length_width_ratio_all,length_axis_all,cell_analysis_all_txt,"cell_area","length_width_ratio","length_axis");


}

else{
    cout<<"Fatal error: no minor mode selected ! (major mode is Experiment)"<<endl;
    exit(-1);
}

}

else if(major_mode=="plot"){

if(minor_mode=="panel"){
    vector<string> plot_files;
    string file_folder = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/";
    string save_panel_name = file_folder + "db_10.00_panel.png";

    //vector<double> parameter_1_plot = {1.00,0.80,0.60,0.40,0.20,0.00};
    vector<double> parameter_1_plot = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    vector<double> parameter_2_plot = {10.00};
    //vector<double >parameter_2_plot = {0.00,0.20,0.40,0.60,0.80,1.00};

    string parameter_1_str_abbr = "y";
    string parameter_2_str_abbr = "db";

    for(int i=0;i<parameter_1_plot.size();i++){
        for(int j=0;j<parameter_2_plot.size();j++){
           if(parameter_1_plot[i]<=parameter_2_plot[j]){
            char plot_name[100];
            sprintf(plot_name, "%s_%.2f_%s_%.2f_k_0_vtk.png",parameter_1_str_abbr.c_str(),parameter_1_plot[i],parameter_2_str_abbr.c_str(),parameter_2_plot[j]);
            string plot_fileName = plot_name;
            plot_fileName = file_folder+plot_fileName;
            plot_files.push_back(plot_fileName);
           }
        }
    }
    gnu_plot::panel_plot_triangle(plot_files,parameter_2_plot.size(),parameter_1_plot.size(),save_panel_name);
}

else if(minor_mode=="panel_with_background"){
    std::vector<string> plot_files;
    std::string file_folder = "/mnt/e/c_file/2023_11_14_TS_new_batch_all/analysis_result/";
    std::string save_panel_name = file_folder + "panel.png";
    std::vector<double> parameter_1_plot = {1.0,0.8,0.6,0.4,0.2,0.0};
    std::vector<double> parameter_2_plot = {0.10,0.50,1,2,3,5,7.5,10};

    std::string parameter_1_str_abbr = "y";
    std::string parameter_2_str_abbr = "db";

    for(int i=0;i<parameter_1_plot.size();i++){
        for(int j=0;j<parameter_2_plot.size();j++){
            char plot_name[100];
            sprintf(plot_name,"%s_%.2f_%s_%.2f_k_0_vtk.png",parameter_1_str_abbr.c_str(),parameter_1_plot[i],parameter_2_str_abbr.c_str(),parameter_2_plot[j]);
            std::string plot_fileName = plot_name;
            plot_fileName = file_folder+plot_fileName;
            plot_files.push_back(plot_fileName);
        }
    }
    gnu_plot::panel_plot(plot_files,parameter_2_plot.size(),parameter_1_plot.size(),save_panel_name);
}

else if(minor_mode=="movie"){
    string image_folder = "/mnt/c/again_and_again/codes/git/plant_vertex_model/y_40_basal_part_similarity/";
    std::string image_name = "contour_basal_cell_number_";
    string file_extension = ".png";
    string output_file_name = image_folder + "contour_basal.avi";

    int frame_rate = 5;
    int start_index= 1000;
    int end_index  = 16000;
    int index_step = 300;

    // Initialize video writer
    cv::Mat frame = cv::imread(image_folder + image_name + to_string(start_index) + file_extension);
    cv::VideoWriter writer;
    int codec = cv::VideoWriter::fourcc('M', 'J', 'P', 'G');  // You can use other codecs according to your needs
    writer.open(output_file_name, codec, frame_rate, frame.size());

    // Check if the writer has been properly initialized
    if (!writer.isOpened()) {
        std::cerr << "Could not open the output video file for writing" << std::endl;
        return -1;
    }

    // Loop over all images and add them to the video
    for (int i = start_index; i <= end_index; i+=index_step) {
        string frame_name = image_folder + image_name + std::to_string(i) + file_extension;
        frame = cv::imread(frame_name);
        cout<<"Loading frame: "<<frame_name<<endl;
        if (frame.empty()) {
            std::cerr << "Failed to load frame " << frame_name << std::endl;
            continue;
        }

        // Write frame to video
        writer.write(frame);
    }
}

else if(minor_mode=="histogram"){
    string data_source = "/mnt/e/UT_confocal_data/23_06_18_TS/2023_06_18_T_S_2-4_real_40x_Cycle/all_normalized_angles.txt";
    string png_output = "/mnt/e/UT_confocal_data/23_06_18_TS/2023_06_18_T_S_2-4_real_40x_Cycle/all_distribution_histogram.png";

    vector<double> EdU_angles = readV::read_vd(data_source);

    int binAngle = 15;
    double minAngle = 0.0;
    double maxAngle = 180.0;
    double binCount = (maxAngle-minAngle)/binAngle;
    vector<int> binAngles(binCount,0);
    for(double angle : EdU_angles){
        int binIndex = static_cast<int>((angle-minAngle)/binAngle);
        ++binAngles[binIndex];
    }
    cout_fout_debug::cout_vector_int(binAngles);
        FILE *pipe = popen("gnuplot -persist", "w");

    if (pipe != NULL) {
        // Set the terminal to png and output to a file
        fprintf(pipe, "set terminal png\n");
        fprintf(pipe, "set output '%s'\n",png_output.c_str());
        fprintf(pipe, "set style fill solid border lc rgb 'black'\n");
        fprintf(pipe, "set yrange [0:20]\n");
        fprintf(pipe, "set xtics '%d' \n",binAngle);

        // Plot with points
        fprintf(pipe, "plot '-' with boxes lw 2 lc rgb 'light-cyan' notitle \n");
        // Output data points
        for(int i = 0; i < binAngles.size(); i++) {
            int binMiddle = minAngle + i*binAngle + binAngle/2;
            fprintf(pipe, "%d %d\n", binMiddle, binAngles[i]);
        }
        
        // End of data
        fprintf(pipe, "e\n");
        
        // Flush the pipe
        fflush(pipe);

        // Close the pipe
        pclose(pipe);
}
}

else if(minor_mode=="Gaussian_diagram"){
    
}

else if(minor_mode=="from_vtk_to_videos"){
    std::string vtkFilePath = "/mnt/e/c_file/2023_10_22_temporal_TS_repeat/y_30_tbias_500/";
    std::string geometrics_record_read_path = vtkFilePath + "geometrics_record.txt"; 
    std::string outputPath = "/mnt/e/c_file/2023_10_22_temporal_TS_repeat/y_30_tbias_500/";
    
    Geometrics_analysis_class ga = readV::read_geo_output_txt(geometrics_record_read_path);
    std::vector<double> simulation_steps = ga.getColumnData("Step");
    for(int i=0;i<(int)simulation_steps.size();i++){
        Organ* p_g = new Organ;
        char str_cell[100], str_line[100],str_png[100];
        sprintf(str_cell,"2dv_cell%.5d00000.vtk",(int)simulation_steps[i]/100000);
        sprintf(str_line,"2dv_line%.5d00000.vtk",(int)simulation_steps[i]/100000);
        sprintf(str_png,"2dv_%.5d.png",(int)simulation_steps[i]/100000);

        cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
        
        string CellVTK = vtkFilePath+str_cell;
        string LineVTK = vtkFilePath+str_line;
        string outputPNG = outputPath+str_png;
        readV::oneVTK(p_g,CellVTK,LineVTK);
        autoVTK::VTKLineCell(LineVTK, CellVTK, outputPNG);
    }
}

else if(minor_mode=="from_vtk_to_cell_number_videos"){
    std::string vtkFilePath = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/y10/0/";
    std::string outputPath = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/y10/0/";
    std::string cell_number_record = "/mnt/e/c_file/2023_10_27_arrest_front_repeat/y10/0/geometrics_record.txt";
    int cell_number_per_step = 100;
    int minimum_cell_number=0;
    int maximum_cell_number=3000;
    
    
    int not_vtk_file_number = 10;

    std::vector<int> cell_number_frame = wangVector::createVector(minimum_cell_number,maximum_cell_number,cell_number_per_step);
    std::vector<int> file_step;
    std::vector<int> stepValues;
    std::vector<int> cellnumberValues;
    std::string line;
    // Open file
    std::ifstream file(cell_number_record);
    if (!file.is_open()) {
        std::cerr << "Unable to open file cell_number_record ("<<cell_number_record<<")." << std::endl;
        return 1;
    }
    // Read file line by line
    while (getline(file, line)) {
        int x, y;

        std::istringstream iss(line);
        if (iss >> x >> y) {  // Try to read the first two doubles from the line
            stepValues.push_back(x);
            cellnumberValues.push_back(y);
            std::cout<<"step "<<x<<"; cell number: "<<y<<std::endl;
        } else {
            std::cerr << "Could not read two doubles from a line" << std::endl;
        }
    }
    file.close();

    for(int i=0;i<cell_number_frame.size();i++){
        int min_error=1000;
        int corresponding_step;
        for(int j=0;j<cellnumberValues.size();j++){
            if(min_error>abs(cell_number_frame[i]-cellnumberValues[j]))
            {
             corresponding_step = stepValues[j];
             min_error=abs(cell_number_frame[i]-cellnumberValues[j]);
            }
        }
        file_step.push_back(corresponding_step);
        std::cout<<"For cell number "<<cell_number_frame[i]<<" corresponding file step size is "<<file_step[i]<<std::endl;
    }


    int length_tmp = vtkFilePath.length();
    char directory_tmp[length_tmp+1];
    strcpy(directory_tmp, vtkFilePath.c_str());
    chdir(directory_tmp);
    int file_index = file_process::getFileNum(vtkFilePath);
    cout<<"file_index: "<<file_index<<endl;
    if(file_index==0){
        cout<<"Error: "<<directory_tmp <<"does not contain any files"<<endl;
        exit(1);
    }
    cout<<"file_index "<<file_index<<endl;
    char str_cell[100], str_line[100];
    int simulation_steps = (file_index-not_vtk_file_number)/2;
    cout<<"simulation_steps "<<simulation_steps<<endl;

    for(int i=0;i<cell_number_frame.size();i++){
        Organ* p_g = new Organ;
        char str_cell[100], str_line[100], str_png[100];
        sprintf(str_cell,"2dv_cell%.5d00000.vtk",file_step[i]/100000);
        sprintf(str_line,"2dv_line%.5d00000.vtk",file_step[i]/100000);
        sprintf(str_png,"cell_number_%d.png",cell_number_frame[i]/100000);
        cout<<"Reading "<<str_cell<<" and "<<str_line<<endl;
        
        string CellVTK = vtkFilePath+str_cell;
        string LineVTK = vtkFilePath+str_line;
        string outputPNG =  outputPath+str_png;
        readV::oneVTK(p_g,CellVTK,LineVTK);
        autoVTK::VTKLineCell_show_cell_number(p_g,0,20);
        //autoVTK::VTKLineCell(LineVTK, CellVTK,outputPNG);
        if(p_g->p_c.size()>maximum_cell_number){
            i=simulation_steps;
            break;
        }
        p_g->~Organ();
    }
}

else{
    cout<<"Fatal error: no minor mode selected ! (major mode is Plot)"<<endl;
    exit(-1);
}
}

else if(major_mode=="test"){
    random_device rnd;
            mt19937 mt(rnd());
            normal_distribution<> normal_axis(90,2);
    for(int i=0;i<100;i++){
        double axisTheta= normal_axis(mt);
        while(axisTheta<-M_PI||axisTheta>2*M_PI){
            axisTheta= normal_axis(mt);
        }
        if(axisTheta<0){
            axisTheta += M_PI;
        }
        if(axisTheta>M_PI){
            axisTheta -= M_PI;
        }
        std::cout<<axisTheta/3.14*180<<std::endl;
    }
            

}

else if(major_mode=="preparation"){
    //generate individual file directory and "batchParameter.txt" by codes
    vector<double> parameter_group_y = {20,25,30,35,40,45,50};
    vector<double> parameter_group_x = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500};
    /*
    double initial_parameter_group_x  = -1.0;
    double parameter_group_x_step = 0.10;
    double terminal_parameter_group_x = 1.0;
    int parameter_group_x_size = 21;
    
    double initial_parameter_group_y  = -1.0;
    double parameter_group_y_step = 0.10;
    double terminal_parameter_group_y = 1.0;
    int parameter_group_y_size = 21;

    for(int i=0; i<parameter_group_x_size; i++){
        double parameter_x_tmp = pow(10.0,initial_parameter_group_x+i*parameter_group_x_step);
        parameter_group_x.push_back(parameter_x_tmp);
    }
    
    
    for(int i=0; i<parameter_group_y_size; i++){
        double parameter_y_tmp = pow(10.0,initial_parameter_group_y+i*parameter_group_y_step);
        parameter_group_y.push_back(parameter_y_tmp);
    }

    
    cout<<endl;
    */

    cout<<"Parameter group x has the following members: ";
    for(int i=0; i<parameter_group_x.size(); i++){
        cout<<parameter_group_x[i]<<"   ";
    }
    cout<<endl;

    cout<<"Parameter group y has the following members: ";
    for(int i=0; i<parameter_group_y.size(); i++){
        cout<<parameter_group_y[i]<<"   ";
    }

    for(int i=0; i<parameter_group_x.size(); i++){
        //char dir_file[100];
        //sprintf(dir_file,"y_%f/",parameter_group_x[i]);
        //mkdir(dir_file);
        stringstream folder_ss, batchParameter_ss;
        folder_ss<<fixed<<setprecision(2)<<"../y_"<<parameter_group_x[i]<<"/";
        mkdir(folder_ss.str().c_str(),0777);
        batchParameter_ss<<folder_ss.str()<<"batchParameter.txt";
        std::ofstream fout(batchParameter_ss.str());
        for(int j=0; j<parameter_group_y.size(); j++){
            if(j<parameter_group_y.size()-1){
                fout<<fixed<<setprecision(4)<<parameter_group_x[i]<<" "<<parameter_group_y[j]<<endl;
            }
            else{
                fout<<fixed<<setprecision(4)<<parameter_group_x[i]<<" "<<parameter_group_y[j];
            }
        }
        fout.close();
    }
    

}

else if(major_mode=="openCV"){
    string image_file = "../opencv_example/6-epi.jpg";
    
    cv::Mat img = cv::imread(image_file);
    //check if the image was loaded successfully
    if(img.empty()){
        std::cerr<<"could not open or find the image." <<endl;
        return -1;
    }

    cv::Mat img_8bit;
    img.convertTo(img_8bit, CV_8U);

    cv::namedWindow("Image", cv::WINDOW_AUTOSIZE);
    cv::imshow("Image", img_8bit);

    cv::waitKey(0); // Wait for a keystroke in the window

}

else{
    cout<<"Fatal error: no major mode selected !"<<endl;
    exit(-1);
}

//record termination time
time_t terminal_time_last = wangSystem::terminal_time();
terminal_time.push_back(terminal_time_last);

output::simulation_log(initial_time,terminal_time);

return 0;
}