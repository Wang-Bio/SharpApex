/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef _DIVISION_H
#define _DIVISION_H

#include "class2dv.h"
#include "vecInoue.h"
#include "wangMath.h"
#include "parameter2dv.h"
#include "geo2dv.h"
#include "force2dv.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <random>
#include <math.h> 


namespace division_frequency{
    void no_control(Organ*);
    void area_control(Organ*);
    void balance_control(Organ*);
    void balance_control_BF_3(Organ*);
    void balance_control_BF_1(Organ*);
    double calcGaussian(double,double,double);
    void Gaussian_control(Organ*,double,double);
    void Gaussian_xy_control(Organ*,double,double,double,double);
    void tag_control(Organ *);
    void temporal_Gaussian_linear_mu(Organ*,double,double);
    void temporal_Gaussian_linear_sigma(Organ*,double,double);
    void temporal_even_to_Gaussian(Organ*);
    void uniform_to_Gaussian_control(Organ*,double,double,double);
    void biregion_frequency_position(Organ*,double,double);
    void biregion_frequency_identity(Organ*,double);
    void biregion_frequency_angles_position(Organ*,double,double);
    void arrest_front(Organ*,double,double,double);
    void arrest_front_biregion(Organ*,double,double,double,double);
    void meristem_position_relative_constant(Organ*,double);
}

namespace division_direction{
    _vec<double> angles(Organ*,int);
    _vec<double> random(Organ*,int);
    _vec<double> epi_anticlinal(Organ*,int);
    _vec<double> epi_periclinal(Organ*,int);
    _vec<double> constant_0(Organ*,int);
    _vec<double> constant_90(Organ*,int);
    _vec<double> yin_distribution(Organ*,int);
    _vec<double> mochizuki_bias(Organ*,int,double,double);
    _vec<double> mochizuki_bias_asymmetrical(Organ*,int,double,double,double,double);
    _vec<double> mochizuki_bias_apical_basal(Organ*,int,double,double,double,double,double);
    _vec<double> temporal_angle_bias(Organ*,int,double,double,double,double);
    _vec<double> Gaussian_bias(Organ*,int,double, double);
    _vec<double> biregion_angles_position(Organ*,int,double,double,double,double,double);
    _vec<double> biregion_angles_gradual(Organ*,int,double,double,double,double,double,double);
    _vec<double> biregion_angles_identity(Organ*,int,double,double,double,double,double);
    _vec<double> biregion_frequency_angles_position(Organ*,int,double,double,double,double,double);
    _vec<double> Gaussian_bias_continuous_Gaussian(Organ*,int,double,double,double,double);
    _vec<double> temporal_biregion_angles(Organ*,int,double,double,double,double,double,double,double);
    _vec<double> temporal_angle_bias_Gaussian(Organ*,int,double,double,double,double,int,int);
}

namespace division{
    void cell_time_initialization(Organ*);
    void cellTimeAdd(Organ *p_g, int increase_index);
    void Global(Organ*);
    void One(Organ*, int);
    void Record(Organ*,int); //have issue in recording frequency_modifier: the new cell's frequency_modifier is set to be 1 by default
    void tag_control(Organ*);

}

#endif