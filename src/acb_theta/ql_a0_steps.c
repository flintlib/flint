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
acb_theta_ql_split(const arb_mat_t cho)
{
    slong g = arb_mat_nrows(cho);
    arb_t cmp;
    slong k;

    arb_init(cmp);

    for (k = g - 1; k >= 1; k--)
    {
        arb_mul_2exp_si(cmp, arb_mat_entry(cho, k - 1, k - 1),
            ACB_THETA_QL_SPLIT);
        if (arb_lt(cmp, arb_mat_entry(cho, k, k)))
        {
            break;
        }
    }

    arb_clear(cmp);
    return k - 1;
}

static void
acb_theta_ql_a0_step(acb_ptr r, acb_srcptr rz, acb_srcptr r0, arb_srcptr dist,
    arb_srcptr dist0, slong k, int has_t, int has_z, slong g, slong prec)
{
    slong n = 1 << g;
    arb_ptr d, d0;
    acb_ptr next;
    slong nb_t = (has_t ? 3 : 1);
    slong nb_r = (has_t ? 2 : 1);
    slong nb_z = (has_z ? 2 : 1);

    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);
    next = _acb_vec_init(nb_t * nb_z * n);

    _arb_vec_scalar_mul_2exp_si(d, dist, n, k);
    _arb_vec_scalar_mul_2exp_si(d0, dist0, n, k);

    if (has_t)
    {
        acb_theta_ql_step_3(next, r, r, r0 + k * nb_r * n, d, d0, g, prec);
        if (has_z)
        {
            acb_theta_ql_step_3(next + nb_t * n, r + nb_t * n, r,
                rz + k * nb_r * n, d, d0, g, prec);
        }
    }
    else
    {
        acb_theta_ql_step_1(next, r, r, r0 + k * nb_r * n, d, d0, g, prec);
        if (has_z)
        {
            acb_theta_ql_step_1(next + nb_t * n, r + nb_t * n, r,
                rz + k * nb_r * n, d, d0, g, prec);
        }
    }
    _acb_vec_set(r, next, nb_t * nb_z * n);

    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
    _acb_vec_clear(next, nb_t * nb_z * n);
}

int
acb_theta_ql_a0_steps(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    arb_srcptr dist0, const acb_mat_t tau, slong guard, slong prec,
    acb_theta_ql_worker_t worker)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    int has_z = !_acb_vec_is_zero(z, g);
    slong nb_t = (has_t ? 3 : 1);
    slong nb_z = (has_z ? 2 : 1);
    acb_mat_t w;
    arb_mat_t cho;
    acb_ptr x, roots;
    arb_ptr new_dist;
    slong d, nb_steps;
    slong k;
    int res = 1;

    acb_mat_init(w, g, g);
    arb_mat_init(cho, g, g);
    x = _acb_vec_init(g);
    new_dist = _arb_vec_init(n);

    acb_theta_eld_cho(cho, tau, ACB_THETA_LOW_PREC);
    d = acb_theta_ql_split(cho);
    nb_steps = acb_theta_ql_nb_steps(cho, d, prec);

    roots = _acb_vec_init(nb_z * nb_t * n * nb_steps);

    /* Get roots */
    res = acb_theta_ql_roots(roots, t, x, dist0, tau, nb_steps, guard, prec);
    if (res && has_z)
    {
        res = acb_theta_ql_roots(roots + nb_t * n * nb_steps, t, z, dist, tau,
            nb_steps, guard, prec);
    }

    if (res)
    {
        /* Call a0_naive at 0 */
        acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist0, n, nb_steps);
        res = acb_theta_ql_a0_split(r, t, x, new_dist, w, d, guard, prec, worker);
    }
    if (res && has_z)
    {
        /* Call a0_naive at z */
        _acb_vec_scalar_mul_2exp_si(x, z, g, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist, n, nb_steps);
        res = acb_theta_ql_a0_split(r + nb_t * n, t, x, new_dist, w, d,
            guard, prec, worker);
    }

    if (res)
    {
        for (k = nb_steps - 1; k >= 0; k++)
        {
            acb_theta_ql_a0_step(r, roots + nb_t * n * nb_steps, roots, dist,
                dist0, k, has_t, has_z, g, prec);
        }
    }

    acb_mat_clear(w);
    arb_mat_clear(cho);
    _acb_vec_clear(x, g);
    _arb_vec_clear(new_dist, n);
    _acb_vec_clear(roots, nb_z * nb_t * n * nb_steps);
    return res;
}
