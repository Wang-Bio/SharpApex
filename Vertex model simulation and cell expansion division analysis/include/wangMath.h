/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef WANGMATH_H
#define WANGMATH_H
#include <stdio.h>
#include <math.h>

#include <vector>
#include <string>
#include <random>
#include "parameter2dv.h"
#include "class2dv.h"

namespace wangMath{
    double mochizuki_bias_sampling_pdf(double,double,double);
    //requires the sampling_pdf
    double mochizuki_bias_single_random_sampling(double,double);

    Quadratic_Solution Quadratic_Equation_Solve(double,double,double);
}

namespace wangVector{
    template <typename T>
    std::vector<T> createVector(T min, T max, T step){
        std::vector<T> result;
        const double epsilon = 1e-10;

        for(T value=min; value<=max+epsilon; value +=step){
            result.push_back(value);
        }
        return result;
    }
}


#endif