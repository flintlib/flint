/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static int
acb_theta_eld_contains_rec(const acb_theta_eld_t E, const slong * pt)
{
    slong d = E->dim;
    slong c = pt[d - 1];
    slong k;

    if (c < (E->min) || c > (E)->max)
    {
        return 0;
    }
    else if (d == 1)
    {
        return 1;
    }
    else if (c >= (E->mid))
    {
        k = c - (E->mid);
        return acb_theta_eld_contains_rec(&E->rchildren[k], pt);
    }
    else
    {
        k = (E->mid) - 1 - c;
        return acb_theta_eld_contains_rec(&E->lchildren[k], pt);
    }
}

int
acb_theta_eld_contains(const acb_theta_eld_t E, const slong * pt)
{
    slong g = E->ambient_dim;
    slong d = E->dim;
    slong k;

    if (acb_theta_eld_nb_pts(E) == 0)
    {
        return 0;
    }

    for (k = d; k < g; k++)
    {
        if (pt[k] != E->last_coords[k - d])
        {
            return 0;
        }
    }

    return acb_theta_eld_contains_rec(E, pt);
}
