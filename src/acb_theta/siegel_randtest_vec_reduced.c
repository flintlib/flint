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
acb_siegel_randtest_vec_reduced(acb_ptr zs, flint_rand_t state, slong nb,
    const acb_mat_t tau, int exact, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_t t;
    slong j;

    /* Sample z in [-1,1] + tau[-1,1] */
    arb_init(t);

    _acb_vec_zero(zs, nb * g);
    for (j = 0; j < nb * g; j++)
    {
        arb_urandom(acb_realref(&zs[j]), state, prec);
        acb_mul_2exp_si(&zs[j], &zs[j], 1);
        acb_sub_si(&zs[j], &zs[j], 1, prec);
    }
    for (j = 0; j < nb; j++)
    {
        acb_mat_vector_mul_col(zs + j * g, tau, zs + j * g, prec);
    }
    for (j = 0; j < nb * g; j++)
    {
        arb_urandom(t, state, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_sub_si(t, t, 1, prec);
        acb_add_arb(&zs[j], &zs[j], t, prec);
    }

    if (exact)
    {
        for (j = 0; j < nb * g; j++)
        {
            acb_get_mid(&zs[j], &zs[j]);
        }
    }

    arb_clear(t);
}
