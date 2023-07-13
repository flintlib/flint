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

static void acb_theta_uql_a0_basecase(acb_ptr th, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong g, slong prec)
{
    acb_mat_t tau_rec;
    slong k, j;

    if (g == 0)
    {
        for (k = 0; k < nb_z; k++)
        {
            acb_one(&th[k]);
        }
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
        acb_theta_ql_a0(th, z, nb_z, tau_rec, prec);        
        acb_mat_clear(tau_rec);
    }
}

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
acb_theta_uql_a0_rec(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau,
    const arb_mat_t cho, const arf_t R2, slong g, slong d, slong prec)
{
    acb_ptr th_rec, z_rec, fac;
    arb_mat_t Yinv;
    slong* ind_z_rec;
    slong nb_z_rec;
    slong n, pt;
    ulong a;
    slong k, j;
    int has_pt;

    flint_printf("(uql_a0_rec) calling with d = %wd, g = %wd, nb_z = %wd\n", d, g, nb_z);

    if (d == g)
    {
        acb_theta_uql_a0_basecase(th, z, nb_z, tau, g, prec);
        return;
    }
    
    /* Recursive call, perhaps doubling nb_z to account for a = 0 and 1 */
    n = 1 << (g - 1);
    th_rec = _acb_vec_init(2 * nb_z * n);
    z_rec = _acb_vec_init(2 * (g - 1) * nb_z);
    nb_z_rec = 0;
    ind_z_rec = flint_malloc(2 * nb_z * sizeof(slong));
    fac = _acb_vec_init(2 * nb_z);
    arb_mat_init(Yinv, g, g);

    /* Set Yinv */
    for (k = 0; k < g; k++)
    {
        for (j = 0; j < g; j++)
        {
            arb_set(arb_mat_entry(Yinv, j, k), acb_imagref(acb_mat_entry(tau, j, k)));
        }
    }
    arb_mat_inv(Yinv, Yinv, prec);

    /* Collect z_rec and cofactors */
    for (k = 0; k < nb_z; k++)
    {        
        for (a = 0; a <= 1; a++)
        {
            /* Get lattice points */
            has_pt = acb_theta_uql_a0_has_pt(&pt, z + k * g, Yinv, cho, R2, a,
                ACB_THETA_ELD_DEFAULT_PREC);
            
            if (has_pt)
            {
                ind_z_rec[2 * k + a] = nb_z_rec;
                
                /* Get new vector z */
                _acb_vec_set(z_rec + (g - 1) * nb_z_rec, z + k * g, g - 1);
                for (j = 0; j < g - 1; j++)
                {
                    acb_addmul_si(z_rec + (g - 1) * nb_z_rec + j,
                        acb_mat_entry(tau, j, g), pt, prec);
                }
                /* Get cofactor */
                acb_mul_si(&fac[2 * k + a], acb_mat_entry(tau, g - 1, g - 1), pt, prec);
                acb_addmul_si(&fac[2 * k + a], &z[k * g + (g - 1)], 2, prec);
                acb_mul_si(&fac[2 * k + a], &fac[2 * k + a], pt, prec);
                acb_exp_pi_i(&fac[2 * k + a], &fac[2 * k + a], prec);
                
                nb_z_rec += 1;
            }
            else
            {
                ind_z_rec[2 * k + a] = -1;
                
                for (j = 0; j < n; j++)
                {
                    acb_zero(&th[k * 2 * n + 2 * j + a]);
                }
            }
        }
    }
    
    /* Recursive call and reconstruct theta values */
    acb_theta_uql_a0_rec(th_rec, z_rec, nb_z_rec, tau, cho, R2, g - 1, d, prec);
    
    flint_printf("th_rec:\n");
    _acb_vec_printd(th_rec, nb_z_rec * n, 10);
    flint_printf("\n");

    for (k = 0; k < nb_z; k++)
    {
        for (a = 0; a <= 1; a++)
        {
            flint_printf("k = %wd, a = %wd, ind = %wd\n", k, a, ind_z_rec[2*k+a]);
            if (ind_z_rec[2 * k + a] == -1)
            {
                continue;
            }
            for (j = 0; j < n; j++)
            {
                acb_mul(&th[k * 2 * n + 2 * j + a], &fac[2 * k + a],
                    &th_rec[ind_z_rec[2 * k + a] * n + j], prec);
            }
        }
    }
    
    _acb_vec_clear(th_rec, 2 * nb_z * n);
    _acb_vec_clear(z_rec, 2 * (g - 1) * nb_z);
    flint_free(ind_z_rec);
    _acb_vec_clear(fac, 2 * nb_z);
    arb_mat_clear(Yinv);
}

void
acb_theta_uql_a0(acb_ptr th, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec)
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
    acb_theta_uql_a0_rec(th, z, nb_z, tau, cho, R2, g, d, prec);

    arb_mat_clear(Y);
    arb_mat_clear(cho);
    arb_clear(pi);
    arf_clear(R2);
    arf_clear(eps);
}
