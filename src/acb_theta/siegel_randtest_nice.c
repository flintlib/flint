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
acb_siegel_randtest_nice(acb_mat_t tau, flint_rand_t state, slong prec)
{
    slong g = arb_mat_nrows(tau);
    acb_t c;
    slong k, j;

    acb_init(c);

    acb_mat_zero(tau);
    for (k = 0; k < g; k++)
    {
        acb_onei(c);
        acb_mul_si(c, c, k + 3, prec);
        acb_mul_2exp_si(c, c, -1);
        acb_set(acb_mat_entry(tau, k, k), c);
    }

    for (k = 0; k < g; k++)
    {
        for (j = k + 1; j < g; j++)
        {
            arb_urandom(acb_realref(c), state, prec);
            arb_sub_si(acb_realref(c), acb_realref(c), 1, prec);
            arb_urandom(acb_imagref(c), state, prec);
            arb_sub_si(acb_imagref(c), acb_imagref(c), 1, prec);
            acb_div_si(c, c, 2*g, prec);
            acb_add(acb_mat_entry(tau, k, j),
                    acb_mat_entry(tau, k, j), c, prec);
            acb_set(acb_mat_entry(tau, j, k), acb_mat_entry(tau, k, j));
        }
    }

    acb_clear(c);
}
