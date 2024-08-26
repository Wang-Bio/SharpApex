#ifndef _GROWTH_H
#define _GROWTH_H

#include "class.h"

extern std::string growth_mode;

extern double homogeneous_growth_rate;
extern double boundary_apical;
extern double boundary_basal;
extern double heterogeneous_growth_rate_x;
extern double heterogeneous_growth_rate_y;
extern double quadratic_mid_growth_rate;


namespace growth{
    void contour(Contour* ct);
    Growth_vector homogeneous_growth(Point_element *pt, Contour *ct);
    Growth_vector biregion_growth(Point_element *pt, Contour *ct);
    Growth_vector biregion_growth_quadratic(Point_element *pt, Contour *ct);
    Growth_vector polarized_growth(Point_element *pt, Contour *ct);
    
}

#endif