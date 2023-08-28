/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_ql_all_sqr(acb_ptr th2, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong guard = ACB_THETA_LOW_PREC;
    int has_z = !_acb_vec_is_zero(z, g);
    slong nb_z = (has_z ? 2 : 1);
    slong nb_t = 1;
    flint_rand_t state;
    acb_mat_t w;
    arb_ptr dist, dist0;
    acb_ptr t, x, th;
    slong j, k;
    int res;

    flint_randinit(state);
    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    dist = _arb_vec_init(n);
    dist0 = _arb_vec_init(n);
    t = _acb_vec_init(g);
    th = _acb_vec_init(n * 3 * nb_z);

    acb_mat_scalar_mul_2exp_si(w, tau, 1);
    _acb_vec_scalar_mul_2exp_si(x, z, g, 1);

    acb_theta_dist_a0(dist, x, w, lp);
    acb_theta_dist_a0(dist0, t, w, lp);

    res = acb_theta_ql_a0(th, t, x, dist0, dist, w, guard, prec);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        nb_t = 3;
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec);
        }
        _acb_vec_scalar_mul_2exp_si(t, t, g, 1);
        res = acb_theta_ql_a0(th, t, x, dist0, dist, w, guard, prec);
        guard += ACB_THETA_LOW_PREC;
    }

    if (!res)
    {
        _acb_vec_indeterminate(th2, n * n);
    }
    else if (has_z)
    {
        acb_theta_ql_dupl(th2, th, th + nb_t * n, dist0, dist, g, prec);
    }
    else
    {
        acb_theta_ql_dupl(th2, th, th, dist0, dist0, g, prec);
    }

    flint_randclear(state);
    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    _arb_vec_clear(dist, n);
    _arb_vec_clear(dist0, n);
    _acb_vec_clear(t, g);
    _acb_vec_clear(th, n * 3 * nb_z);
}
