/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_agm_sqrt_entry(acb_t r, const acb_t a, const acb_t root, slong prec)
{
    acb_t y1, y2;

    acb_init(y1);
    acb_init(y2);

    acb_sqrts(y1, y2, a, prec);

    if (acb_overlaps(root, y1) && acb_overlaps(root, y2))
    {
        acb_indeterminate(r);
    }
    else if (acb_overlaps(root, y1))
    {
        acb_set(r, y1);
    }
    else if (acb_overlaps(root, y2))
    {
        acb_set(r, y2);
    }
    else
    {
        flint_printf("(agm_sqrt) Error: no overlap\n");
        acb_printd(a, 10); flint_printf("\n");
        acb_printd(root, 10); flint_printf("\n");
        flint_abort();
        acb_indeterminate(r);
    }

    acb_clear(y1);
    acb_clear(y2);
}

void
acb_theta_agm_sqrt(acb_ptr r, acb_srcptr a, acb_srcptr roots, slong nb, slong prec)
{
    slong k;

    for (k = 0; k < nb; k++)
    {
        acb_theta_agm_sqrt_entry(&r[k], &a[k], &roots[k], prec);
    }
}
