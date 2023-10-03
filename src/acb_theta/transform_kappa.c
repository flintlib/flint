/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_siegel_sqrtdet_i(acb_t r, const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    slong prec = ACB_THETA_LOW_PREC;
    acb_mat_t tau;
    acb_t x;
    fmpz_t re, im;
    int done = 0;

    acb_mat_init(tau, g, g);
    acb_mat_onei(tau);
    acb_init(x);
    fmpz_init(re);
    fmpz_init(im);

    while (!done)
    {
        prec *= 2;
        acb_siegel_cocycle_det(x, mat, tau, prec);
        done = arb_get_unique_fmpz(re, acb_realref(r))
            && arb_get_unique_fmpz(im, acb_imagref(r));
    }
    arb_set_fmpz(acb_realref(r), re);
    arb_set_fmpz(acb_imagref(r), im);
    acb_sqrts(r, x, x, prec);

    acb_mat_clear(tau);
    fmpz_clear(re);
    fmpz_clear(im);
}

static slong
get_power_of_zeta8(const acb_t x)
{
    acb_t y;
    arb_t abs;
    slong prec = ACB_THETA_LOW_PREC;
    slong k;

    acb_init(y);
    arb_init(abs);

    for (k = 0; k < 8; k++)
    {
        acb_one(y);
        acb_mul_2exp_si(y, y, -2);
        acb_mul_si(y, y, -k, prec);
        acb_exp_pi_i(y, y, prec);
        acb_mul(y, y, x, prec);
        arb_abs(abs, acb_imagref(y));
        if (arb_lt(abs, acb_realref(y)))
        {
            /* y is in correct quadrant */
            break;
        }
    }
    acb_sub_si(y, y, 1, prec);

    if (k < 8 && !acb_contains_zero(y))
    {
        flint_printf("(acb_theta_transform_k) Error: not a power of zeta8\n");
        flint_printf("k = %wd, y:\n", k);
        acb_printd(y, 10);
        flint_printf("\n");
        flint_abort();
    }
    else if (k == 8)
    {
        k = -1; /* insufficient precision */
    }

    acb_clear(y);
    arb_clear(abs);
    return k;
}

slong
acb_theta_transform_kappa(const fmpz_mat_t mat)
{
    slong g = sp2gz_dim(mat);
    fmpz_mat_t inv;
    acb_mat_t tau;
    acb_ptr z;
    acb_t scal1, scal2, t;
    fmpz_t eps;
    ulong ab;
    slong kappa = -1;
    slong prec = ACB_THETA_LOW_PREC;

    fmpz_mat_init(inv, 2 * g, 2 * g);
    acb_mat_init(tau, g, g);
    z = _acb_vec_init(g);
    fmpz_init(eps);
    acb_init(scal1);
    acb_init(scal2);
    acb_init(t);

    sp2gz_inv(inv, mat);
    ab = acb_theta_transform_char(eps, inv, 0);
    acb_theta_transform_char(eps, mat, ab);

    while (kappa == -1)
    {
        acb_mat_onei(tau);
        acb_theta_naive_00(scal1, z, 1, tau, prec);

        acb_siegel_sqrtdet_i(t, mat);
        acb_siegel_transform(tau, mat, tau, prec);
        acb_theta_naive_fixed_ab(scal2, ab, z, 1, tau, prec);

        acb_mul(scal1, scal1, t, prec);
        acb_set_fmpz(t, eps);
        acb_mul_2exp_si(t, t, -2);
        acb_exp_pi_i(t, t, prec);
        acb_mul(scal1, scal1, t, prec);
        acb_div(scal1, scal2, scal1, prec);

        kappa = get_power_of_zeta8(scal1);
        prec *= 2;
    }

    fmpz_mat_clear(inv);
    acb_mat_clear(tau);
    _acb_vec_clear(z, g);
    fmpz_clear(eps);
    acb_clear(scal1);
    acb_clear(scal2);
    acb_clear(t);
    return kappa;
}
