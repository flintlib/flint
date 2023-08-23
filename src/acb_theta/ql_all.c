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
    acb_ptr roots, new_z, th_a0;
    arb_ptr new_dist0, new_dist;
    slong hprec;
    slong k;
    int res = 1;

    acb_mat_init(new_tau, g, g);
    roots = _acb_vec_init(n * n);
    new_z = _acb_vec_init(g);
    th_a0 = _acb_vec_init(n * nb_z * nb_t);
    new_dist0 = _arb_vec_init(n);
    new_dist = _arb_vec_init(n);

    /* Collect roots: we only need theta_{a,b}(z + t, tau) */
    _acb_vec_add(new_z, z, t, g, prec);
    for (a = 0; a < n; a++)
    {
        hprec = guard + acb_theta_dist_addprec(&dist[a]);
        acb_theta_naive_fixed_a(roots + a * n, a << g, new_z, 1, tau, hprec);

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
        acb_theta_duplication(th, th_a0, th_a0 + (nb_z * nb_t - 1) * n,
            new_dist0, new_dist, g, prec);
        acb_theta_agm_sqrt(th, th, roots, n * n, prec);

        if (has_t)
        {
            /* Get theta_{a,b}(z, tau) from division */
            acb_theta_duplication(aux, th_a0 + n, th_a0 + (3 * nb_z - 2) * n,
                new_dist0, new_dist, g, prec);
            for (k = 0; k < n * n; k++)
            {
                acb_div(&th[k], &aux[k], &th[k], prec);
            }
        }
    }

    acb_mat_clear(new_tau, g, g);
    _acb_vec_clear(roots, n * n);
    _acb_vec_clear(new_z, g);
    _acb_vec_clear(th_a0, n * nb_z * nb_t);
    _arb_vec_clear(new_dist0, n);
    _arb_vec_clear(new_dist, n);
    return res;
}

void
acb_theta_ql_all(acb_ptr th, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong guard = ACB_THETA_LOW_PREC;
    slong nb_t = 1;
    flint_rand_t state;
    arb_ptr dist, dist0;
    acb_ptr t;
    slong j;
    int res;

    flint_randinit(state);
    dist = _arb_vec_init(n);
    dist0 = _arb_vec_init(n);
    t = _acb_vec_init(g);

    acb_theta_dist_a0(dist, z, tau, lp);
    acb_theta_dist_a0(dist0, t, tau, lp);

    res = acb_theta_ql_all_with_t(th, t, z, dist0, dist, tau, guard, prec);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        nb_t = 3;
        res = acb_theta_ql_a0_with_t(th, t, z, dist0, dist, tau, guard, prec);
        guard += ACB_THETA_LOW_PREC;
    }
    if (!res)
    {
        _acb_vec_indeterminate(th, n * n);
    }

    flint_randclear(state);
    _arb_vec_clear(dist, n);
    _arb_vec_clear(dist0, n);
    _acb_vec_clear(t, g);
}
