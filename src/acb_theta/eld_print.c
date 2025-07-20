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
acb_theta_eld_print(const acb_theta_eld_t E)
{
    slong d = E->dim;
    slong g = E->ambient_dim;
    slong k;

    for (k = 0; k < g - d; k++)
    {
        flint_printf("  ");
    }
    flint_printf("Slice (...");
    for (k = 0; k < g - d; k++)
    {
        flint_printf(", %wd", E->last_coords[k]);
    }
    flint_printf("): from %wd to %wd (mid: %wd)\n", E->min, E->max, E->mid);
    if (d > 1)
    {
        for (k = 0; k < (E->nr); k++)
        {
            acb_theta_eld_print(&E->rchildren[k]);
        }
        for (k = 0; k < (E->nl); k++)
        {
            acb_theta_eld_print(&E->lchildren[k]);
        }
    }
}
