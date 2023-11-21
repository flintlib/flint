/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = acb_mat_nrows(tau);
    slong s = n_randint(state, g + 1);
    slong n = n_randint(state, FLINT_MAX(1, mag_bits));
    fmpz_mat_t mat;
    arb_t test;
    int r = 0;
    slong j, k;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    arb_init(test);

    for (k = 0; (k < 10) && !r; k++)
    {
        acb_siegel_randtest(tau, state, prec, 2);
        acb_siegel_reduce(mat, tau, prec);
        acb_siegel_transform(tau, mat, tau, prec);
        r = acb_siegel_is_reduced(tau, -1, prec);
    }
    if (!r)
    {
        acb_mat_onei(tau);
    }

    for (j = s; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            arb_mul_2exp_si(acb_imagref(acb_mat_entry(tau, j, k)),
                acb_imagref(acb_mat_entry(tau, j, k)), n);
            arb_mul_2exp_si(acb_imagref(acb_mat_entry(tau, k, j)),
                acb_imagref(acb_mat_entry(tau, k, j)), n);
        }
    }

    fmpz_mat_clear(mat);
    arb_clear(test);
}
