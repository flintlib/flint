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
acb_theta_eld_border(slong * pts, const acb_theta_eld_t E)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong k, i;

    if (d == 1)
    {
        pts[0] = (E->min) - 1;
        pts[g] = (E->max) + 1;
        for (k = 1; k < g; k++)
        {
            pts[k] = E->last_coords[k - d];
            pts[k + g] = E->last_coords[k - d];
        }
    }
    else /* d > 1 */
    {
        i = 0;
        for (k = 0; k < (E->nr); k++)
        {
            acb_theta_eld_border(&pts[i], &E->rchildren[k]);
            i += g * acb_theta_eld_nb_border(&E->rchildren[k]);
        }
        for (k = 0; k < (E->nl); k++)
        {
            acb_theta_eld_border(&pts[i], &E->lchildren[k]);
            i += g * acb_theta_eld_nb_border(&E->lchildren[k]);
        }
    }
}
