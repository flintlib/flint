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
acb_theta_eld_init(acb_theta_eld_t E, slong d, slong g)
{
    FLINT_ASSERT(d >= 1 && d <= g);

    acb_theta_eld_dim(E) = d;
    acb_theta_eld_ambient_dim(E) = g;
    E->last_coords = flint_malloc((g - d) * sizeof(slong));
    E->rchildren = NULL;
    acb_theta_eld_nr(E) = 0;
    E->lchildren = NULL;
    acb_theta_eld_nl(E) = 0;
    E->box = flint_malloc(d * sizeof(slong));
}
