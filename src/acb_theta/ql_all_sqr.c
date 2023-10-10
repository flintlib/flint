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
acb_theta_ql_all_sqr_red(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong guard = ACB_THETA_LOW_PREC;
    int hasz = !_acb_vec_is_zero(z, g);
    slong nbz = (hasz ? 2 : 1);
    slong nbt = 1;
    flint_rand_t state;
    acb_mat_t w;
    arb_ptr d, d0;
    acb_ptr t, x, th;
    slong j, k;
    int res;

    flint_randinit(state);
    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);
    t = _acb_vec_init(g);
    th = _acb_vec_init(n * 3 * nbz);

    acb_mat_scalar_mul_2exp_si(w, tau, 1);
    _acb_vec_scalar_mul_2exp_si(x, z, g, 1);

    acb_theta_dist_a0(d, x, w, lp);
    acb_theta_dist_a0(d0, t, w, lp);

    res = acb_theta_ql_a0(th, t, x, d0, d, w, guard, prec);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        nbt = 3;
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec);
        }
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        res = acb_theta_ql_a0(th, t, x, d0, d, w, guard, prec);
        guard += ACB_THETA_LOW_PREC;
    }

    if (!res)
    {
        _acb_vec_indeterminate(th2, n * n);
    }
    else if (hasz)
    {
        acb_theta_ql_dupl(th2, th, th + nbt * n, d0, d, g, prec);
    }
    else
    {
        acb_theta_ql_dupl(th2, th, th, d0, d0, g, prec);
    }

    flint_randclear(state);
    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
    _acb_vec_clear(t, g);
    _acb_vec_clear(th, n * 3 * nbz);
}

void
acb_theta_ql_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    acb_mat_t w;
    acb_ptr new_z, aux;
    acb_t c;
    arb_t u, v;
    arf_t b;
    slong s, j, k;
    ulong ab, a0, a1, b0, b1, fixed_a1;

    acb_init(c);
    arb_init(u);
    arb_init(v);
    arf_init(b);
    new_z = _acb_vec_init(g);

    s = acb_theta_ql_reduce(new_z, c, u, &fixed_a1, z, tau, prec);
    acb_sqr(c, c, prec);
    arb_sqr(u, u, prec);

    if (s == -1)
    {
        _acb_vec_zero(th2, n2);
        for (ab = 0; ab < n2; ab++)
        {
            acb_add_error_arb(&th2[ab], u);
        }
    }
    else
    {
        acb_mat_init(w, s, s);
        aux = _acb_vec_init(1 << (2 * s));

        for (j = 0; j < s; j++)
        {
            for (k = 0; k < s; k++)
            {
                acb_set(acb_mat_entry(w, j, k), acb_mat_entry(tau, j, k));
            }
        }

        if (acb_is_finite(c))
        {
            if (s > 0)
            {
                acb_theta_ql_all_sqr_red(aux, new_z, w, prec);
            }
            else
            {
                acb_one(&aux[0]);
            }
            _acb_vec_scalar_mul(aux, aux, 1 << (2 * s), c, prec);
        }
        else
        {
            _acb_vec_indeterminate(aux, 1 << (2 * s));
        }

        for (ab = 0; ab < n2; ab++)
        {
            /* Write ab as a0 a1 b0 b1 */
            a0 = ab >> (g + (g - s));
            a1 = (ab >> g) % (1 << (g - s));
            b0 = (ab >> (g - s)) % (1 << s);
            b1 = ab % (1 << (g - s));

            if (a1 == fixed_a1)
            {
                acb_set(&th2[ab], &aux[(a0 << s) + b0]);
                if (acb_theta_char_dot(a1, b1, g - s) % 2 == 1)
                {
                    acb_neg(&th2[ab], &th2[ab]);
                }
                acb_abs(v, &th2[ab], prec);
                arb_mul(v, v, u, prec);
                arb_get_ubound_arf(b, v, prec);
                arb_set_arf(v, b);
                arb_sqrt(v, v, prec);
                arb_mul_2exp_si(v, v, 1);
                acb_add_error_arb(&th2[ab], v);
            }
            else
            {
                acb_zero(&th2[ab]);
                acb_add_error_arb(&th2[ab], u);
            }
        }

        acb_mat_clear(w);
        _acb_vec_clear(aux, 1 << (2 * s));
    }

    _acb_vec_clear(new_z, g);
    acb_clear(c);
    arb_clear(u);
    arb_clear(v);
    arf_clear(b);
}