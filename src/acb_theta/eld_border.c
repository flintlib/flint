/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_eld_border(slong * pts, const acb_theta_eld_t E)
{
    slong d = acb_theta_eld_dim(E);
    slong g = acb_theta_eld_ambient_dim(E);
    slong nr = acb_theta_eld_nr(E);
    slong nl = acb_theta_eld_nl(E);
    slong max = acb_theta_eld_max(E);
    slong min = acb_theta_eld_min(E);
    slong k, i;

    if (d == 1)
    {
        pts[0] = min - 1;
        pts[g] = max + 1;
        for (k = 1; k < g; k++)
        {
            pts[k] = acb_theta_eld_coord(E, k);
            pts[k + g] = acb_theta_eld_coord(E, k);
        }
    }
    else /* d > 1 */
    {
        i = 0;
        for (k = 0; k < nr; k++)
        {
            acb_theta_eld_border(&pts[i], acb_theta_eld_rchild(E, k));
            i += g * acb_theta_eld_nb_border(acb_theta_eld_rchild(E, k));
        }
        for (k = 0; k < nl; k++)
        {
            acb_theta_eld_border(&pts[i], acb_theta_eld_lchild(E, k));
            i += g * acb_theta_eld_nb_border(acb_theta_eld_lchild(E, k));
        }
    }
}
