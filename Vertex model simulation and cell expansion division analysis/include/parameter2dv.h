/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef PARAMETER2DV_H
#define PARAMETER2DV_H

#include "class2dv.h"
#include "IO2dv.h"
#include "vecInoue.h"
#include "vtk2dv.h"

#include <string>
#include <fstream>
#include <iostream>

using namespace std;

//debug_mode
constexpr const char *debug_force_record ="OFF"; //when "ON", force will be recorded for each T_vtkoutput

//potential_energy_mode 
constexpr const char *potential_energy_mode = "simple"; //1. simple, suggested by Mochizuki-sensei, mostly used; 2. L_std

//list of parameters
extern int repeat_time;
extern int checkend_tmp;

//parameters about force
extern double sigma_L;
extern double sigma_O;
extern double kappa_S;
extern double S_std;
extern double eta;
extern double L_std;

//parameters about time and cell division
extern double standard_cell_period_length;
extern double F_apparent;
extern double F_modifier; //modify cell_division_threshold
extern int T_division_check;
extern int T_vtkoutput;
extern int T_cell_time;
extern double delta_time ; //originally set as 5e-3
extern int step_end ;
extern int end_cell_number ; //end when cell number reaches this number 

//"area_control"
extern double area_control_lower_limit;
extern double area_control_slope;

//"balance_control"
extern double BF;

//"Gaussian_control"
extern double gau_sigma;
extern double gau_mu;

extern double temporal_Gaussian_a;
extern double temporal_Gaussian_b;

//"temporal_Gaussian_linear_mu" parameters
extern double temporal_Gaussian_initial_mu;
extern double temporal_Gaussian_terminal_mu;

//"temporal_Gaussian_linear_sigma" parameters
extern double temporal_Gaussian_initial_sigma;
extern double temporal_Gaussian_terminal_sigma;

// "Gaussian_xy" parameters
extern double gau_mu_x;
extern double gau_sigma_x;
extern double gau_mu_y;
extern double gau_sigma_y;

//"uniform_to_Gaussian"
extern double uniform_to_Gaussian_transition_time;

extern int temporal_t1;
extern int temporal_t2;

//"mochizuki_bias" parameters
extern double mochizuki_bias_phi;
extern double mochizuki_bias_beta;

//"mochizuki_bias_asymmetrical" parameters
extern double mochizuki_bias_beta_left;
extern double mochizuki_bias_beta_right;
extern double mochizuki_bias_phi_left;
extern double mochizuki_bias_phi_right;

//"mochizuki_bias_apical_basal" parameters
extern double mochizuki_bias_beta_apical;
extern double mochizuki_bias_beta_basal;
extern double mochizuki_bias_phi_apical;
extern double mochizuki_bias_phi_basal;
extern double angle_bias_y_boundary;

//"temporal_angle_bias" parameters
extern double temporal_bias_beta_initial;
extern double temporal_bias_beta_terminal;
extern double temporal_bias_phi_initial;
extern double temporal_bias_phi_terminal;

//"Gaussian_bias" parameters
extern double Gaussian_bias_phi;
extern double Gaussian_bias_beta;

//"biregion_angles_position" parameters
extern double biregion_angles_position_apical_gaussian_beta;
extern double biregion_angles_position_basal_gaussian_beta;
extern double biregion_angles_position_apical_gaussian_phi;
extern double biregion_angles_position_basal_gaussian_phi;
extern double biregion_angles_position_y_boundary;

//"biregion_angles_identity" parameters
extern double biregion_angles_identity_apical_gaussian_beta;
extern double biregion_angles_identity_basal_gaussian_beta;
extern double biregion_angles_identity_apical_gaussian_phi;
extern double biregion_angles_identity_basal_gaussian_phi;
extern double biregion_angles_identity_y_boundary;

//"biregion_frequency_position" parameters
extern double biregion_frequency_position_y_boundary;
extern double biregion_frequency_position_relative_frequency;

//"biregion_frequency_identity" parameters
extern double biregion_identity_y_boundary;
extern double biregion_identity_relative_frequency;

//"biregion_frequency_angles_position" parameters
extern double biregion_frequency_angles_position_apical_gaussian_beta;
extern double biregion_frequency_angles_position_basal_gaussian_beta;
extern double biregion_frequency_angles_position_apical_gaussian_phi;
extern double biregion_frequency_angles_position_basal_gaussian_phi;
extern double biregion_frequency_angles_position_relative_frequency;
extern double biregion_frequency_angles_position_y_boundary;

//"Gaussian_bias_continuous_Gaussian" parameters
extern double Gaussian_bias_beta_A;
extern double Gaussian_bias_beta_mu;
extern double Gaussian_bias_beta_sigma;
//extern double Gaussian_bias_phi;

//"temporal_biregion_angles" parameters
extern double y_boundary_initial;
extern double y_boundary_change;
extern double y_boundary_terminal;
extern double temporal_biregion_angles_apical_gaussian_beta;
extern double temporal_biregion_angles_basal_gaussian_beta;
extern double temporal_biregion_angles_apical_gaussian_phi;
extern double temporal_biregion_angles_basal_gaussian_phi;

//"arrest_front" parameters
extern double y_arrest_front;
extern double t_arrest_front;
extern double k_arrest_front;

//"simga_O_spatial" parameters
extern double sigma_O_spatial_max;
extern double sigma_O_spatial_min;
extern double sigma_O_spatial_y_boundary;

//"S_std_spatial" parameters
extern double S_std_base;
extern double S_std_apex;

//"temporal_angle_bias_Gaussian"
extern double bias_sigma_initial;
extern double bias_sigma_terminal;
extern double bias_mu_initial;
extern double bias_mu_terminal;
extern int bias_t_init;
extern int bias_t_term;

//meristem_position_relative_constant
extern double meristem_position_relative_constant_y;

//biregion_angles_gradual
extern double biregion_angles_gradual_apical_gaussian_beta;
extern double biregion_angles_gradual_basal_gaussian_beta;
extern double biregion_angles_gradual_apical_gaussian_phi;
extern double biregion_angles_gradual_basal_gaussian_phi;
extern double biregion_angles_gradual_y1_boundary; //boundary between apical and joint
extern double biregion_angles_gradual_y2_boundary; //boundary between joint and basal

//arrest_front_biregion
extern double arrest_front_biregion_y;
extern double arrest_front_biregion_x;
extern double arrest_front_biregion_t;
extern double arrest_front_biregion_k;

//list of files for input and output
extern string parameterInputFile;
extern string parameterRecordFile;
extern string initialOrganFile;
extern string modeFile;

//list of modes used
extern string major_mode;
extern string minor_mode;

extern string division_control;
extern string in_division_direction;
extern string epi_division_direction;
extern string temporal_mode;

extern string mechanics_mode;


namespace parameter{
    void read_mode(string);
    void read(string);
    void record();
    void batchRead(parameterList*);
}

namespace termination{
    bool checkEnd(Organ*);
}

#endif