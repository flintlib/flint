/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_randtest_compact(acb_mat_t tau, flint_rand_t state, int exact, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong bits = n_randint(state, 4);
    arb_t y;
    slong j, k;
    int res;

    arb_init(y);

    res = 0;
    while(!res)
    {
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        arb_sub_si(y, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 100, prec);
        res = arb_is_negative(y);
    }

    if (exact)
    {
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_get_mid(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k));
                acb_set(acb_mat_entry(tau, k, j), acb_mat_entry(tau, j, k));
            }
        }
    }

    arb_clear(y);
}
