/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static slong
acb_theta_uql_cut(const arb_mat_t cho, const arf_t R2, slong prec)
{
    slong g = arb_mat_nrows(cho);
    arb_t x;
    arf_t bound;
    slong res = g;

    arb_init(x);
    arf_init(bound);

    while (res > 0)
    {
        arb_sqr(x, arb_mat_entry(cho, res - 1, res - 1), prec);
        arb_mul_2exp_si(x, x, -2);
        arb_get_ubound_arf(bound, x, prec);
        if (arf_cmp(bound, R2) < 0)
        {
            break;
        }
        else
        {
            res -= 1;
        }
    }

    arb_clear(x);
    arf_clear(bound);
    return res;
}

static void acb_theta_uql_a0_basecase(acb_ptr th, acb_srcptr z,
    const acb_mat_t tau, slong g, slong prec)
{
    acb_mat_t tau_rec;
    slong k, j;

    if (g == 0)
    {
        acb_one(&th[0]);
    }
    else
    {
        acb_mat_init(tau_rec, g, g);
        for (k = 0; k < g; k++)
        {
            for (j = 0; j < g; j++)
            {
                acb_set(acb_mat_entry(tau_rec, j, k), acb_mat_entry(tau, j, k));
            }
        }
        acb_theta_ql_a0(th, z, 1, tau_rec, prec);
        acb_mat_clear(tau_rec);
    }
}

/* If return 0, leave pt unchanged */
static int
acb_theta_uql_a0_has_pt(slong *pt, acb_srcptr z, const arb_mat_t Yinv,
    const arb_mat_t cho, const arf_t R2, ulong a, slong prec)
{
    slong g = arb_mat_nrows(Yinv);
    arb_ptr v;
    arb_t c, x;
    arf_t rad;
    slong min, mid, max;
    int res;

    v = _arb_vec_init(g);
    arb_init(c);
    arb_init(x);
    arf_init(rad);

    /* Lattice is Z^g with offset Y^-1 z and Gram matrix Y */
    /* Get center */
    _acb_vec_get_imag(v, z, g);
    arb_mat_vector_mul_col(v, Yinv, v, prec);
    arb_neg(c, &v[g - 1]);
    
    /* Get radius */
    arb_set_arf(x, R2);
    arb_sqrt(x, x, prec);
    arb_div(x, x, arb_mat_entry(cho, g - 1, g - 1), prec);
    arb_get_ubound_arf(rad, x, prec);

    /* Get points */
    acb_theta_eld_interval(&min, &mid, &max, c, rad, a, prec);
    if (min < max)
    {
        /* This should not happen. */
        flint_printf("(uql_a0) Error: found several lattice points\n");
        flint_abort(); /* Replace by indeterminate */
    }
    else if (min == max)
    {
        res = 1;
        *pt = min;
    }
    else
    {
        res = 0;
    }

    _arb_vec_clear(v, g);
    arb_clear(c);
    arb_clear(x);
    arf_clear(rad);
    return res;
}

static void
acb_theta_uql_a0_rec(acb_ptr th, acb_srcptr z, const acb_mat_t tau,
    const arb_mat_t cho, const arf_t R2, slong g, slong d, slong prec)
{
    acb_ptr th_rec, z_rec;
    acb_t fac;
    arb_mat_t Yinv;
    slong n, pt;
    slong k, j;
    int has_pt_0, has_pt_1;
    ulong a;

    flint_printf("(uql_a0_rec) calling with d = %wd, g = %wd\n", d, g);

    if (d == g)
    {
        acb_theta_uql_a0_basecase(th, z, tau, g, prec);
        return;
    }
    
    /* Recursive call */
    n = 1 << (g - 1);
    th_rec = _acb_vec_init(n);
    z_rec = _acb_vec_init(g - 1);
    arb_mat_init(Yinv, g, g);
    acb_init(fac);

    /* Set Yinv */
    for (k = 0; k < g; k++)
    {
        for (j = 0; j < g; j++)
        {
            arb_set(arb_mat_entry(Yinv, j, k), acb_imagref(acb_mat_entry(tau, j, k)));
        }
    }
    arb_mat_inv(Yinv, Yinv, prec);

    has_pt_0 = acb_theta_uql_a0_has_pt(&pt, z, Yinv, cho, R2, 0,
        ACB_THETA_ELD_DEFAULT_PREC);
    has_pt_1 = acb_theta_uql_a0_has_pt(&pt, z, Yinv, cho, R2, 1,
        ACB_THETA_ELD_DEFAULT_PREC);
    
    if (has_pt_0 && has_pt_1)
    {
        flint_printf("(uql_a0) Error: two points_n");
        flint_abort();
    }
    else if (!has_pt_0 && !has_pt_1)
    {
        _acb_vec_zero(th, 2 * n);
    }
    else /* One point */
    {
        a = has_pt_0 ? 0 : 1;
        _acb_vec_set(z_rec, z, g - 1);
        for (j = 0; j < g - 1; j++)
        {
            acb_addmul_si(&z_rec[j], acb_mat_entry(tau, j, g), pt, prec);
        }
        acb_mul_si(fac, acb_mat_entry(tau, g - 1, g - 1), pt, prec);
        acb_addmul_si(fac, &z[g - 1], 2, prec);
        acb_mul_si(fac, fac, pt, prec);
        acb_exp_pi_i(fac, fac, prec);

        acb_theta_uql_a0_rec(th_rec, z_rec, tau, cho, R2, g - 1, d, prec);
        
        flint_printf("th_rec:\n");
        _acb_vec_printd(th_rec, n, 10);
        flint_printf("\n");

        for (j = 0; j < n; j++)
        {
            acb_mul(&th[2 * j + a], &th_rec[j], fac, prec);
            acb_zero(&th[2 * j + ((a + 1) % 2)]);
        }
    }
    
    _acb_vec_clear(th_rec, n);
    _acb_vec_clear(z_rec, g - 1);
    arb_mat_clear(Yinv);
    acb_clear(fac);
}

void
acb_theta_uql_a0(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t Y;
    arb_mat_t cho;
    arb_t pi;
    arf_t R2;
    arf_t eps;
    slong ord = 0;
    slong d;

    arb_mat_init(Y, g, g);
    arb_mat_init(cho, g, g);
    arb_init(pi);
    arf_init(R2);
    arf_init(eps);

    acb_mat_get_imag(Y, tau);
    arb_const_pi(pi, prec);
    arb_mat_scalar_mul_arb(cho, Y, pi, prec);
    arb_mat_cho(cho, cho, prec);
    arb_mat_transpose(cho, cho);

    arf_one(eps);
    arf_mul_2exp_si(eps, eps, -prec);
    acb_theta_naive_radius(R2, Y, ord, eps, prec);

    d = acb_theta_uql_cut(cho, R2, ACB_THETA_ELD_DEFAULT_PREC);
    acb_theta_uql_a0_rec(th, z, tau, cho, R2, g, d, prec);

    arb_mat_clear(Y);
    arb_mat_clear(cho);
    arb_clear(pi);
    arf_clear(R2);
    arf_clear(eps);
}
