/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

static void
acb_theta_agm_sqrt_entry(acb_t res, const acb_t a, const acb_t rt, slong prec)
{
    acb_t y1, y2;
    int t1, t2;

    acb_init(y1);
    acb_init(y2);

    acb_sqrts(y1, y2, a, prec);
    t1 = acb_overlaps(rt, y1);
    t2 = acb_overlaps(rt, y2);

    if (t1 && t2)
    {
        acb_union(res, y1, y2, prec);
    }
    else if (t1)
    {
        acb_set(res, y1);
    }
    else if (t2)
    {
        acb_set(res, y2);
    }
    else
    {
        acb_indeterminate(res);
    }

    acb_clear(y1);
    acb_clear(y2);
}

void
acb_theta_agm_sqrt(acb_ptr res, acb_srcptr a, acb_srcptr rts, slong nb, slong prec)
{
    slong k;

    for (k = 0; k < nb; k++)
    {
        acb_theta_agm_sqrt_entry(&res[k], &a[k], &rts[k], prec);
    }
}
