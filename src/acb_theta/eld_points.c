/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_eld_points(slong * pts, const acb_theta_eld_t E)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong k, j, i;

    if (d == 1)
    {
        i = 0;
        for (k = (E->min); k <= (E->max); k++)
        {
            pts[i] = k;
            for (j = 1; j < g; j++)
            {
                pts[i + j] = E->last_coords[j - d];
            }
            i += g;
        }
    }
    else /* d > 1 */
    {
        i = 0;
        for (k = 0; k < (E->nr); k++)
        {
            acb_theta_eld_points(&pts[i], &E->rchildren[k]);
            i += g * ((&E->rchildren[k])->nb_pts);
        }
        for (k = 0; k < (E->nl); k++)
        {
            acb_theta_eld_points(&pts[i], &E->lchildren[k]);
            i += g * ((&E->lchildren[k])->nb_pts);
        }
    }
}
