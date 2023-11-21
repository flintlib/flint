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

#define ACB_THETA_QL_TRY 100

static void
acb_theta_ql_dupl(acb_ptr th2, acb_srcptr th0, acb_srcptr th,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr v;
    ulong a, b;

    v = _acb_vec_init(n);

    for (a = 0; a < n; a++)
    {
        _acb_vec_set(v, th, n);
        for (b = 0; b < n; b++)
        {
            if (acb_theta_char_dot(a, b, g) % 2 == 1)
            {
                acb_neg(&v[b], &v[b]);
            }
        }
        acb_theta_agm_mul_tight(v, th0, v, d0, d, g, prec);
        for (b = 0; b < n; b++)
        {
            acb_set(&th2[b * n + a], &v[b]);
        }
    }
    _acb_vec_scalar_mul_2exp_si(th2, th2, n * n, g);

    _acb_vec_clear(v, n);
}


static int
acb_theta_ql_all_with_t(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hasz = !_acb_vec_is_zero(z, g);
    int hast = !_acb_vec_is_zero(t, g);
    slong nbz = (hasz ? 2 : 1);
    slong nbt = (hast ? 3 : 1);
    acb_mat_t new_tau;
    acb_ptr rts, new_z, th_a0, aux;
    arb_ptr new_d0, new_d;
    slong hprec;
    slong k, a;
    int res = 1;

    acb_mat_init(new_tau, g, g);
    rts = _acb_vec_init(n * n);
    new_z = _acb_vec_init(g);
    th_a0 = _acb_vec_init(n * nbz * nbt);
    aux = _acb_vec_init(n * n);
    new_d0 = _arb_vec_init(n);
    new_d = _arb_vec_init(n);

    /* Collect roots: we only need theta_{a,b}(z + t, tau) */
    _acb_vec_add(new_z, z, t, g, prec);

    for (a = 0; (a < n) && res; a++)
    {
        hprec = guard + acb_theta_dist_addprec(&d[a]);
        acb_theta_naive_fixed_a(rts + a * n, a, new_z, 1, tau, hprec);
        for (k = 0; (k < n) && res; k++)
        {
            /* Ignore theta constants if z = t = 0 */
            if (acb_contains_zero(&rts[a * n + k])
                && (hasz || hast || acb_theta_char_is_even(a * n + k, g)))
            {
                res = 0;
            }
        }
    }

    /* Get ql_a0 at 2z, t, 2tau */
    if (res)
    {
        acb_mat_scalar_mul_2exp_si(new_tau, tau, 1);
        _acb_vec_scalar_mul_2exp_si(new_z, z, g, 1);
        _arb_vec_scalar_mul_2exp_si(new_d, d, n, 1);
        _arb_vec_scalar_mul_2exp_si(new_d0, d0, n, 1);

        res = acb_theta_ql_a0(th_a0, t, new_z, new_d0, new_d, new_tau, guard, prec);
    }

    if (res)
    {
        /* Get theta_{a,b}(z + t, tau) from square roots */
        acb_theta_ql_dupl(th, th_a0, th_a0 + (nbz * nbt - 1) * n,
            new_d0, new_d, g, prec);
        acb_theta_agm_sqrt(th, th, rts, n * n, prec);

        if (hast)
        {
            /* Get theta_{a,b}(z, tau) from division */
            acb_theta_ql_dupl(aux, th_a0 + n, th_a0 + (3 * nbz - 2) * n,
                new_d0, new_d, g, prec);
            for (k = 0; k < n * n; k++)
            {
                acb_div(&th[k], &aux[k], &th[k], prec);
            }
        }
    }

    /* Set odd theta constants to zero */
    if (!hasz)
    {
        for (a = 0; a < n * n; a++)
        {
            if (!acb_theta_char_is_even(a, g))
            {
                acb_zero(&th[a]);
            }
        }
    }

    acb_mat_clear(new_tau);
    _acb_vec_clear(rts, n * n);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(th_a0, n * nbz * nbt);
    _acb_vec_clear(aux, n * n);
    _arb_vec_clear(new_d0, n);
    _arb_vec_clear(new_d, n);
    return res;
}

static void
acb_theta_ql_all_red(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong guard = ACB_THETA_LOW_PREC;
    slong nb_der = acb_theta_jet_nb(2, g);
    flint_rand_t state;
    arb_ptr d, d0;
    acb_mat_t tau_mid;
    acb_ptr t, z_mid, dth;
    arb_t err;
    arf_t e;
    slong j, k;
    int hasz = !_acb_vec_is_zero(z, g);
    int res;

    flint_randinit(state);
    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);
    acb_mat_init(tau_mid, g, g);
    t = _acb_vec_init(g);
    z_mid = _acb_vec_init(g);
    dth = _acb_vec_init(n * n * nb_der);
    arb_init(err);
    arf_init(e);

    acb_theta_dist_a0(d, z, tau, lp);
    acb_theta_dist_a0(d0, t, tau, lp);

    /* Get midpoints; ql_all_with_t is expected to lose guard + g bits of precision */
    arf_one(e);
    arf_mul_2exp_si(e, e, -prec - guard - g);
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_get_mid(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau, j, k));
            acb_add_error_arf(acb_mat_entry(tau_mid, j, k), e);
            acb_set(acb_mat_entry(tau_mid, k, j), acb_mat_entry(tau_mid, j, k));
        }
        acb_get_mid(&z_mid[j], &z[j]);
        if (hasz)
        {
            acb_add_error_arf(&z_mid[j], e);
        }
    }

    res = acb_theta_ql_all_with_t(th, t, z_mid, d0, d, tau_mid,
        guard, prec + guard + g);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec + guard + g);
        }
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        res = acb_theta_ql_all_with_t(th, t, z_mid, d0, d, tau_mid,
            guard, prec + guard + g);
        guard += ACB_THETA_LOW_PREC;
    }

    if (!res)
    {
        _acb_vec_indeterminate(th, n * n);
    }
    else
    {
        acb_theta_jet_naive_all(dth, z, tau, 2, ACB_THETA_LOW_PREC);
        for (j = 0; j < n * n; j++)
        {
            acb_theta_jet_error_bounds(err, z, tau, dth + j * nb_der, 0, ACB_THETA_LOW_PREC);
            acb_add_error_arb(&th[j], err);
        }
    }

    flint_randclear(state);
    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(t, g);
    _acb_vec_clear(z_mid, g);
    _acb_vec_clear(dth, n * n * nb_der);
    arb_clear(err);
    arf_clear(e);
}

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
acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    acb_mat_t tau0;
    acb_ptr new_z, aux;
    acb_t c;
    arb_t u, v;
    arf_t b;
    slong s;
    slong * n1;
    ulong ab, a0, a1, b0, b1, fixed_a1;

    acb_init(c);
    arb_init(u);
    arb_init(v);
    arf_init(b);
    new_z = _acb_vec_init(g);
    n1 = flint_malloc(g * sizeof(slong));

    s = acb_theta_ql_reduce(new_z, c, u, n1, z, tau, prec);
    if (sqr)
    {
        acb_sqr(c, c, prec);
        arb_sqr(u, u, prec);
    }

    if (s == -1)
    {
        _acb_vec_zero(th, n2);
        for (ab = 0; ab < n2; ab++)
        {
            acb_add_error_arb(&th[ab], u);
        }
    }
    else
    {
        fixed_a1 = acb_theta_char_get_a(n1, g - s);
        acb_mat_window_init(tau0, tau, 0, 0, s, s);
        aux = _acb_vec_init(1 << (2 * s));

        if (acb_is_finite(c))
        {
            if (s > 0 && !sqr)
            {
                acb_theta_ql_all_red(aux, new_z, tau0, prec);
            }
            else if (s > 0)
            {
                acb_theta_ql_all_sqr_red(aux, new_z, tau0, prec);
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

            if (a1 != fixed_a1)
            {
                acb_zero(&th[ab]);
            }
            else
            {
                acb_mul_i_pow_si(&th[ab], &aux[(a0 << s) + b0],
                    (sqr ? 2 : 1) * acb_theta_char_dot_slong(b1, n1, g - s));
                if (sqr)
                {
                    acb_abs(v, &th[ab], prec);
                    arb_mul(v, v, u, prec);
                    arb_get_ubound_arf(b, v, prec);
                    arb_set_arf(v, b);
                    arb_sqrt(v, v, prec);
                    arb_mul_2exp_si(v, v, 1);
                    acb_add_error_arb(&th[ab], v);
                }
            }
            acb_add_error_arb(&th[ab], u);
        }

        acb_mat_window_clear(tau0);
        _acb_vec_clear(aux, 1 << (2 * s));
    }

    _acb_vec_clear(new_z, g);
    acb_clear(c);
    arb_clear(u);
    arb_clear(v);
    arf_clear(b);
    flint_free(n1);
}
