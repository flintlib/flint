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
acb_siegel_randtest_fund(acb_mat_t tau, flint_rand_t state, slong prec)
{
    slong g = arb_mat_nrows(tau);
    arf_t rad;
    acb_t err;
    acb_t c;
    slong k, j;

    arf_init(rad);
    acb_init(err);
    acb_init(c);

    arf_one(rad);
    arf_div_si(rad, rad, 2 * g, prec, ARF_RND_FLOOR);

    acb_mat_zero(tau);
    for (k = 0; k < g; k++)
    {
        acb_onei(c);
        acb_mul_si(c, c, k + 3, prec);
        acb_mul_2exp_si(c, c, -1);
        acb_set(acb_mat_entry(tau, k, k), c);
    }

    acb_zero(c);
    for (k = 0; k < g; k++)
    {
        for (j = k + 1; j < g; j++)
        {
            acb_randtest_disk(err, c, rad, state, prec);
            acb_add(acb_mat_entry(tau, k, j),
                    acb_mat_entry(tau, k, j), err, prec);
            acb_set(acb_mat_entry(tau, j, k), acb_mat_entry(tau, k, j));
        }
    }

    arf_clear(rad);
    acb_clear(err);
    acb_clear(c);
}
