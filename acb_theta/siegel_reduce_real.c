/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_siegel_reduce_real(fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong j, k;
    fmpz_t c;

    if (!acb_mat_is_finite(tau))
    {
        fmpz_mat_zero(mat);
        return;
    }

    fmpz_init(c);

    fmpz_mat_one(mat);
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            arf_get_fmpz(c, arb_midref(acb_realref(acb_mat_entry(tau, j, k))),
                         ARF_RND_NEAR);
            fmpz_neg(fmpz_mat_entry(mat, j, k + g), c);
        }
        for (k = 0; k < j; k++)
        {
            fmpz_set(fmpz_mat_entry(mat, j, k + g),
                     fmpz_mat_entry(mat, k, j + g));
        }
    }

    fmpz_clear(c);
}
