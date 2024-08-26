#ifndef INITIAL2DV_H
#define INITIAL2DV_H

#include "Organ2dv.h"
#include "class2dv.h"
#include "force2dv.h"
#include "geo2dv.h"
#include "division.h"
#include "IO2dv.h"

extern string real_organ_contour_imagej_txt;

namespace initialization{
    void organ(Organ*);
    void organ_continue(Organ*);
}


#endif