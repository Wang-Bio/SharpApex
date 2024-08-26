/*************************************************************************************************/
// 2D Vertex Model for Plant Morphogenesis
// Original Author: Yasuhiro Inoue (inoue.yasuhiro.4n@kyoto-u.ac.jp), for animal morphogenesis simulation
// Modified by: Zining Wang (wangzining16@mails.ucas.ac.cn), for plant morphogenesis simulation based on cell division patterns
// Reference: Kinoshita, A., Naito, M., Wang, Z., Inoue, Y., Mochizuki, A., & Tsukaya, H. (2022). Position of meristems and the angles of the cell division plane regulate the uniqueness of lateral organ shape. Development, 149(23), dev199773.
/*********************************************************************************************/

#ifndef FORCE2DV_H
#define FORCE2DV_H

#include "class2dv.h"
#include "vecInoue.h"
#include "parameter2dv.h"
#include "geo2dv.h"
#include <iostream>
#include <cstring>

namespace force {
    void line_elastic_force(Organ*);
    void cell_elastic_force(Organ*);
    void calcForceMotion(Organ*);
    void forceShapeInitiation(Organ*,int);

    //for specific mechanics mode
    //1. S_std_spatial
    void S_std_spatial_geo(Organ*);
}

#endif