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
ql_use_step_1(acb_ptr r, slong nb_z, acb_srcptr roots, arb_srcptr dist,
    slong g, slong prec)
{
    slong n = 1 << g;
    slong k;

    for (k = 1; k < nb_z; k++)
    {
        acb_theta_ql_step(r + k * n, r + k * n, r, roots + k * n,
            dist + k * n, g, prec);
    }
    acb_theta_ql_step(r, r, r, roots, dist, g, prec);
}

static void
ql_use_step_3(acb_ptr r, slong nb_z, acb_srcptr roots, arb_srcptr dist,
    slong g, slong prec)
{
    slong n = 1 << g;
    slong k;
    acb_ptr next;

    next = _acb_vec_init(3 * nb_z * n);

    for (k = 0; k < nb_z; k++)
    {
        acb_theta_ql_step_aux(next + 3 * k * n, r + 3 * k * n, r,
            roots + 2 * k * n, dist + k * n, g, prec);
    }

    _acb_vec_set(r, next, 3 * nb_z * n);
    _acb_vec_clear(next, 3 * nb_z * n);
}

int
acb_theta_ql_use_steps(acb_ptr r, acb_srcptr t, acb_srcptr z, slong nb_z,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec,
    acb_theta_ql_worker_t worker_d)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    arb_mat_t cho;
    acb_ptr x, roots;
    arb_ptr new_dist;    
    slong d, nb_steps;
    slong k;
    int res = 1;

    acb_mat_init(w, g, g);
    arb_mat_init(cho, g, g);
    x = _acb_vec_init(nb_z * g);
    roots = _acb_vec_init(3 * nb_z * n);
    new_dist = _arb_vec_init(nb_z * n);
    
    acb_theta_eld_cho(cho, tau, prec);
    d = acb_theta_ql_cut(cho);
    nb_steps = acb_theta_ql_new_nb_steps(cho, d, prec);

    /* Get roots */
    res = acb_theta_ql_any_roots(roots, t, z, nb_z, dist, tau, nb_steps, guard, prec);

    if (res)
    {
        /* Call use_naive */
        acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
        _acb_vec_scalar_mul_2exp_si(x, z, nb_z * g, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_dist, dist, nb_z * n, nb_steps);
        for (k = 0; (k < nb_z) && res; k++)
        {
            res = acb_theta_ql_use_naive(r, t, x + k * g, new_dist + k * n, w,
                d, prec, worker_d);
        }
    }

    if (res)
    {        
        /* Make steps */
        for (k = nb_steps - 1; k >= 0; k++)
        {
            if (_acb_vec_is_zero(t, g))
            {
                /* This is a problem because roots aren't organized like this */
                ql_use_step_1(r, nb_z, roots + k * n * nb_z, dist, g, prec);
            }
            else
            {
                ql_use_step_3(r, nb_z, roots + 2 * k * n * nb_z, dist, g, prec);
            }
            _arb_vec_scalar_mul_2exp_si(new_dist, new_dist, n, -1);
        }
    }

    acb_mat_clear(w);
    arb_mat_clear(cho);
    _acb_vec_clear(x, nb_z * g);
    _acb_vec_clear(roots, 3 * nb_z * n);
    _arb_vec_clear(new_dist, nb_z * n);
    return res;
}
