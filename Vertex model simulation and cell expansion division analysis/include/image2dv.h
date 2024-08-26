#ifndef _IMAGEDV_H
#define _IMAGE2DV_H

#include "class2dv.h"
#include <opencv2/opencv.hpp>

class Point2dv;
class Contour2dv;
class Polygon2dv;
class Image2dv;

class Point2dv{
 public:
    int x;
    int y;
    int pixel_value=255;
    bool checked=0;
    int polygon_index=-1;

    Point2dv(int x, int y, int pixel_value) : x(x), y(y), pixel_value(pixel_value) {} 
};

class Contour2dv{
 public:
    std::vector<Point2dv> pt;
};

class Polygon2dv{
 public:
    std::vector<Point2dv> pt;
    double area;
    double length;
    double width;
    double angles;
};

class Image2dv{
 public:
    std::vector<std::vector<Point2dv>> pt;
    std::vector<Polygon2dv> pg;
    int height;
    int width;
};

namespace Image_2dv{
    void read_image(std::string, Image2dv*);
    void draw_image_from_points(Image2dv*);

    Point2dv start_point_in_polygon(Image2dv*);
    void add_white_pixel_to_polygon(Image2dv*,int,int,int);
    bool check_contour_finding_end(Image2dv*);

    void find_polygon_from_image(Image2dv*);
    void find_contours_from_polygon(Image2dv*);
    void statistics_of_polygons(Image2dv*);

    void removeSmallerObjects(cv::Mat&,cv::Mat&);
}

#endif