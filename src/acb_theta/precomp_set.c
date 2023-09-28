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
acb_theta_precomp_set(acb_theta_precomp_t D, acb_srcptr z,
    const acb_mat_t tau, const acb_theta_eld_t E, slong prec)
{
    slong g = acb_theta_eld_ambient_dim(E);
    arb_t pi;
    acb_t c, dc, ddc;
    slong k, j;
    slong nb_pow;

    if (acb_theta_eld_nb_pts(E) == 0)
    {
        return;
    }

    arb_init(pi);
    acb_init(c);
    acb_init(dc);
    acb_init(ddc);

    arb_const_pi(pi, prec);

    /* Set matrix of exponentials */
    /* Matrix has exp(i pi (1 + delta_{j,k}) tau_{j,k}) in upper triangle */
    for (k = 0; k < g; k++)
    {
        for (j = k; j < g; j++)
        {
            acb_mul_arb(c, acb_mat_entry(tau, k, j), pi, prec);
            acb_mul_onei(c, c);
            if (k != j)
            {
                acb_mul_2exp_si(c, c, 1);
            }
            acb_exp(c, c, prec);
            acb_set(acb_mat_entry(acb_theta_precomp_exp_mat(D), k, j), c);
        }
    }

    /* Set indices */
    D->indices[0] = 0;
    for (k = 0; k < g; k++)
    {
        nb_pow = acb_theta_eld_box(E, k) + 1;
        D->indices[k + 1] = D->indices[k] + nb_pow;
    }

    /* Init and set square powers; addition chains unnecessary */
    /* Contain exp(i pi j^2 tau_{k,k}) for j up to box */
    D->sqr_powers = _acb_vec_init(D->indices[g]);
    for (k = 0; k < g; k++)
    {
        acb_one(c);
        acb_set(dc, acb_mat_entry(acb_theta_precomp_exp_mat(D), k, k));
        acb_sqr(ddc, dc, prec);
        for (j = 0; j <= acb_theta_eld_box(E, k); j++)
        {
            acb_set(acb_theta_precomp_sqr_pow(D, k, j), c);
            acb_mul(c, c, dc, prec);
            acb_mul(dc, dc, ddc, prec);
        }
    }

    /* Set exponentials of z */
    /* Contain exp(2 i pi z_j) */
    for (k = 0; k < acb_theta_precomp_nb_z(D); k++)
    {
        for (j = 0; j < g; j++)
        {
            acb_mul_2exp_si(acb_theta_precomp_exp_z(D, k, j), &z[k * g + j], 1);
            acb_exp_pi_i(acb_theta_precomp_exp_z(D, k, j),
                acb_theta_precomp_exp_z(D, k, j), prec);
        }
    }

    arb_clear(pi);
    acb_clear(c);
    acb_clear(dc);
    acb_clear(ddc);
}
