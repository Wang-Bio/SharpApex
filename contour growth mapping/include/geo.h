#ifndef _GEO_H
#define _GEO_H

#include "class.h"

namespace geo{
    void contour_center(Contour *ct);
    void contour_min_max(Contour *ct);
    void continuity_check(Contour* ct);

    void normalized_contour(Contour* ct);
    void sampling_contour(Contour* ct);

}

#endif
