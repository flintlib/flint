/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "acb_theta.h"

acb_theta_ctx_z_struct *
acb_theta_ctx_z_vec_init(slong nb, slong g)
{
    acb_theta_ctx_z_struct * res;
    slong k;

    FLINT_ASSERT(nb >= 0);
    FLINT_ASSERT(g >= 1);

    res = flint_malloc(nb * sizeof(acb_theta_ctx_z_struct));
    for (k = 0; k < nb; k++)
    {
        acb_theta_ctx_z_init(&res[k], g);
    }

    return res;
}
