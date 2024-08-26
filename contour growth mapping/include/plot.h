#ifndef _PLOT_H
#define _PLOT_H


#include "class.h"
#include <opencv2/opencv.hpp>

namespace plot{
    void save_contour_txt(std::string txt_file_name, Contour *ct);
    void save_contour_png(std::string png_file_name, Contour *ct);
    void save_contour_png(std::string png_file_name, Contour *ct,double range_x1, double range_x2, double range_y1, double range_y2);

    void input_real_contour(std::string real_contour_file, Contour* ct);

    void create_image_panel(std::vector<std::string> images, int nRows, int nCols, std::string);
}


#endif