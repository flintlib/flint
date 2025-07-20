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

void
acb_theta_ctx_z_vec_clear(acb_theta_ctx_z_struct * vec, slong nb)
{
    slong k;

    FLINT_ASSERT(nb >= 0);
    for (k = 0; k < nb; k++)
    {
        acb_theta_ctx_z_clear(&vec[k]);
    }
    flint_free(vec);
}
