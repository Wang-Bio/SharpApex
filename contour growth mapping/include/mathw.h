#ifndef _MATHW_H
#define _MATHW_H

#include "class.h"
#include <cmath>

double f_upper(double x,double k){
    return std::sin(x+M_PI/2.0)-k*x;
}

double df_upper(double x, double k){
    return std::cos(x+M_PI/2.0)-k;
}

double newtonRaphson_upper(double initialGuess, double k, int maxIter = 1000, double tolerance = 1e-7) {
    double x = initialGuess;
    for (int i = 0; i < maxIter; ++i) {
        double fx = f_upper(x, k);
        double dfx = df_upper(x, k);

        // Avoid division by zero
        if (std::abs(dfx) < tolerance) {
            return x;
        }

        double x_new = x - fx / dfx;
        
        // Check for convergence
        if (std::abs(x_new - x) < tolerance) {
            return x_new;
        }
        x = x_new;
    }
    return x; // Return the last computed value
}

double f_lower(double x,double k){
    return std::sin(x-M_PI/2.0)-k*x;
}

double df_lower(double x, double k){
    return std::cos(x-M_PI/2.0)-k;
}

double newtonRaphson_lower(double initialGuess, double k, int maxIter = 1000, double tolerance = 1e-7) {
    double x = initialGuess;
    for (int i = 0; i < maxIter; ++i) {
        double fx = f_lower(x, k);
        double dfx = df_lower(x, k);

        // Avoid division by zero
        if (std::abs(dfx) < tolerance) {
            return x;
        }

        double x_new = x - fx / dfx;
        
        // Check for convergence
        if (std::abs(x_new - x) < tolerance) {
            return x_new;
        }
        x = x_new;
    }
    return x; // Return the last computed value
}

#endif