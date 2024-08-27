/*
    Copyright (C) 2023 Jean Kieffer

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
acb_siegel_randtest_vec_reduced(acb_ptr z, flint_rand_t state,
    const acb_mat_t tau, int exact, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t t;
    slong j;

    /* Sample z in [-2,2] + tau[-2,2] */
    arb_init(t);

    _acb_vec_zero(z, g);
    for (j = 0; j < g; j++)
    {
        arb_urandom(acb_realref(&z[j]), state, prec);
        acb_mul_2exp_si(&z[j], &z[j], 2);
        acb_sub_si(&z[j], &z[j], 2, prec);
    }
    acb_mat_vector_mul_col(z, tau, z, prec);
    for (j = 0; j < g; j++)
    {
        arb_urandom(t, state, prec);
        arb_mul_2exp_si(t, t, 2);
        arb_sub_si(t, t, 2, prec);
        acb_add_arb(&z[j], &z[j], t, prec);
    }

    if (exact)
    {
        for (j = 0; j < g; j++)
        {
            acb_get_mid(&z[j], &z[j]);
        }
    }

    arb_clear(t);
}
