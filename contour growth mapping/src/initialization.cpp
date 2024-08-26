#include "../include/initialization.h"

namespace initial{
    void circle(Contour* ct,int element_number,double radius,double center_x,double center_y){
        //add points by polar coordinates
        ct->initial_radius=radius;
        for(int i=0;i<element_number;i++){
            Point_element* pt_tmp = new Point_element;
            double angles_tmp = i*6.2831852/element_number;
            pt_tmp->x = center_x + radius*std::cos(angles_tmp);
            pt_tmp->y = center_y + radius*std::sin(angles_tmp);
            ct->pt.push_back(pt_tmp);
        }
    }
}