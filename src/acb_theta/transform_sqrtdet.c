/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* Use theta relation at low precision */

static void
acb_theta_transform_sqrtdet_lowprec(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong prec = ACB_THETA_LOW_PREC;
    slong kappa = acb_theta_transform_kappa(mat);
    slong b = -1;
    acb_mat_t w;
    acb_ptr z, th;
    acb_t t;
    fmpz_t eps;
    slong k;
    ulong ab;

    acb_mat_init(w, g, g);
    z = _acb_vec_init(g);
    th = _acb_vec_init(n);
    acb_init(t);
    fmpz_init(eps);

    while (b < 0)
    {
        /* Find b such that theta_{0,b}(gamma tau) is nonzero */
        prec *= 2;
        acb_mat_get_mid(w, tau);
        acb_siegel_transform(w, mat, w, prec);
        acb_theta_naive_0b(th, z, 1, w, prec);
        for (k = 0; k < n; k++)
        {
            if (!acb_contains_zero(&th[k]))
            {
                b = k;
                break;
            }
        }
    }

    ab = acb_theta_transform_char(eps, mat, b);
    acb_zero(r);

    while (acb_contains_zero(r))
    {
        acb_theta_naive_fixed_ab(th, b, z, 1, w, prec);
        acb_theta_naive_fixed_ab(r, ab, z, 1, tau, prec);
        acb_set_fmpz(t, eps);
        acb_add_si(t, t, kappa, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_exp_pi_i(t, t, prec);
        acb_mul(r, r, t, prec);
        acb_div(r, &th[b], r, prec);
        prec *= 2;
    }

    acb_mat_clear(w);
    _acb_vec_clear(z, g);
    _acb_vec_clear(th, n);
    acb_clear(t);
    fmpz_clear(eps);
}

void
acb_theta_transform_sqrtdet(acb_t r, const fmpz_mat_t mat, const acb_mat_t tau, slong prec)
{
    acb_t x, y1, y2;

    acb_init(x);
    acb_init(y1);
    acb_init(y2);

    acb_siegel_cocycle_det(r, mat, tau, prec);
    acb_theta_transform_sqrtdet_lowprec(x, mat, tau);
    acb_sqrts(y1, y2, r, prec);

    if (acb_overlaps(y1, x) && acb_overlaps(y2, x))
    {
        acb_union(r, y1, y2, prec);
    }
    else if (acb_overlaps(y1, x))
    {
        acb_set(r, y1);
    }
    else if (acb_overlaps(y2, x))
    {
        acb_set(r, y2);
    }
    else
    {
        flint_printf("(acb_theta_transform_sqrtdet) Error: no overlap\n");
        flint_abort();
    }

    acb_clear(x);
    acb_clear(y1);
    acb_clear(y2);
}
