#ifndef _GNUPLOT_H
#define _GNUPLOT_H

#include "class2dv.h"
#include "Vertex2dv.h"


namespace gnuplot{
    void vv(const std::vector<Vertex> &);
    void vv(const std::vector<Vertex> &,const std::string &);
    void vv(const std::vector<double> &, const std::vector<double> &, const std::string &);
    void curvature_plot(vector<double>);

    void boxplot(const std::vector<double>&, const std::string &);
}

#endif