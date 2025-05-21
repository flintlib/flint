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
acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)
{
    FLINT_ASSERT(d >= 1 && d <= g);

    E->dim = d;
    E->ambient_dim = g;
    E->last_coords = flint_malloc((g - d) * sizeof(slong));
    E->rchildren = NULL;
    E->nr = 0;
    E->lchildren = NULL;
    E->nl = 0;
    E->box = flint_malloc(d * sizeof(slong));
}
