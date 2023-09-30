/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* todo: move out? */
static int
_acb_vec_contains_zero(acb_srcptr v, slong n)
{
    slong k;

    for (k = 0; k < n; k++)
    {
        if (acb_contains_zero(&v[k]))
        {
            return 1;
        }
    }

    return 0;
}

int acb_theta_ql_all_use_naive(slong g, slong prec)
{
    if (g <= 2)
    {
        return (prec <= 1500);
    }
    else
    {
        return 0;
    }
}

static int
acb_theta_ql_all_with_t(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_z = !_acb_vec_is_zero(z, g);
    int has_t = !_acb_vec_is_zero(t, g);
    slong nb_z = (has_z ? 2 : 1);
    slong nb_t = (has_t ? 3 : 1);
    acb_mat_t new_tau;
    acb_ptr roots, new_z, th_a0, aux;
    arb_ptr new_dist0, new_dist;
    slong hprec;
    slong k, a;
    int res = 1;

    acb_mat_init(new_tau, g, g);
    roots = _acb_vec_init(n * n);
    new_z = _acb_vec_init(g);
    th_a0 = _acb_vec_init(n * nb_z * nb_t);
    aux = _acb_vec_init(n * n);
    new_dist0 = _arb_vec_init(n);
    new_dist = _arb_vec_init(n);

    /* Collect roots: we only need theta_{a,b}(z + t, tau) */
    _acb_vec_add(new_z, z, t, g, prec);

    for (a = 0; a < n; a++)
    {
        hprec = guard + acb_theta_dist_addprec(&dist[a]);
        acb_theta_naive_fixed_a(roots + a * n, a, new_z, 1, tau, hprec);

        if (_acb_vec_contains_zero(roots + a * n, n))
        {
            res = 0;
            break;
        }
    }

    /* Get ql_a0 at 2z, t, 2tau */
    if (res)
    {
        acb_mat_scalar_mul_2exp_si(new_tau, tau, 1);
        _acb_vec_scalar_mul_2exp_si(new_z, z, g, 1);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist, n, 1);
        _arb_vec_scalar_mul_2exp_si(new_dist0, dist0, n, 1);

        res = acb_theta_ql_a0(th_a0, t, new_z, new_dist0, new_dist, new_tau, guard, prec);
    }

    if (res)
    {
        /* Get theta_{a,b}(z + t, tau) from square roots */
        acb_theta_ql_dupl(th, th_a0, th_a0 + (nb_z * nb_t - 1) * n,
            new_dist0, new_dist, g, prec);
        acb_theta_agm_sqrt(th, th, roots, n * n, prec);

        if (has_t)
        {
            /* Get theta_{a,b}(z, tau) from division */
            acb_theta_ql_dupl(aux, th_a0 + n, th_a0 + (3 * nb_z - 2) * n,
                new_dist0, new_dist, g, prec);
            for (k = 0; k < n * n; k++)
            {
                acb_div(&th[k], &aux[k], &th[k], prec);
            }
        }
    }

    acb_mat_clear(new_tau);
    _acb_vec_clear(roots, n * n);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(th_a0, n * nb_z * nb_t);
    _acb_vec_clear(aux, n * n);
    _arb_vec_clear(new_dist0, n);
    _arb_vec_clear(new_dist, n);
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
    arb_ptr dist, dist0;
    acb_mat_t tau_mid;
    acb_ptr t, z_mid, dth;
    arb_t err;
    arf_t e;
    slong j, k;
    int has_z = !_acb_vec_is_zero(z, g);
    int res;

    flint_randinit(state);
    dist = _arb_vec_init(n);
    dist0 = _arb_vec_init(n);
    acb_mat_init(tau_mid, g, g);
    t = _acb_vec_init(g);
    z_mid = _acb_vec_init(g);
    dth = _acb_vec_init(n * n * nb_der);
    arb_init(err);
    arf_init(e);

    acb_theta_dist_a0(dist, z, tau, lp);
    acb_theta_dist_a0(dist0, t, tau, lp);

    /* Get midpoints; ql_all_with_t is expected to lose guard+g bits of precision */
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
        if (has_z)
        {
            acb_add_error_arf(&z_mid[j], e);
        }
    }

    res = acb_theta_ql_all_with_t(th, t, z_mid, dist0, dist, tau_mid,
        guard, prec + guard + g);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec + guard + g);
        }
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        res = acb_theta_ql_all_with_t(th, t, z_mid, dist0, dist, tau_mid,
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
    _arb_vec_clear(dist, n);
    _arb_vec_clear(dist0, n);
    acb_mat_clear(tau_mid);
    _acb_vec_clear(t, g);
    _acb_vec_clear(z_mid, g);
    _acb_vec_clear(dth, n * n * nb_der);
    arb_clear(err);
    arf_clear(e);
}

void
acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    acb_mat_t w;
    acb_ptr x, aux;
    acb_t c;
    arb_t u;
    slong d, j, k;
    ulong ab, a0, a1, b0;

    acb_init(c);
    arb_init(u);
    x = _acb_vec_init(g);

    d = acb_theta_ql_reduce(x, c, u, z, tau, prec);

    acb_mat_init(w, d, d);
    aux = _acb_vec_init(1 << (2 * d));

    for (j = 0; j < d; j++)
    {
        for (k = 0; k < d; k++)
        {
            acb_set(acb_mat_entry(w, j, k), acb_mat_entry(tau, j, k));
        }
    }

    if (acb_is_finite(c))
    {
        if (d > 0 && acb_theta_ql_all_use_naive(g, prec))
        {
            acb_theta_naive_all(aux, x, 1, w, prec);
        }
        else if (d > 0)
        {
            acb_theta_ql_all_red(aux, x, w, prec);
        }
        else
        {
            acb_one(&aux[0]);
        }
        _acb_vec_scalar_mul(aux, aux, 1 << (2 * d), c, prec);
    }
    else
    {
        _acb_vec_indeterminate(aux, 1 << (2 * d));
    }

    for (ab = 0; ab < n2; ab++)
    {
        /* Write ab as a0 a1 b0 b1 */
        a0 = ab >> (g + (g - d));
        a1 = (ab >> g) % (1 << (g - d));
        b0 = (ab >> (g - d)) % (1 << d);

        if (a1 == 0)
        {
            acb_set(&th[ab], &aux[(a0 << d) + b0]);
        }
        else
        {
            acb_zero(&th[ab]);
        }
        acb_add_error_arb(&th[ab], u);
    }

    _acb_vec_clear(x, g);
    acb_clear(c);
    arb_clear(u);
    acb_mat_clear(w);
    _acb_vec_clear(aux, 1 << (2 * d));
}
