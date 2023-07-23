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
acb_theta_ql_all_sqr(acb_ptr r, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong lp = ACB_THETA_LOW_PREC;
    slong guard = ACB_THETA_LOW_PREC;
    flint_rand_t state;
    arb_mat_t cho;
    arb_ptr dist;
    acb_ptr t;
    slong k, j;
    slong nb_z = 1; /* Adjust if z is zero */
    int res = 0;

    flint_randinit(state);
    arb_mat_init(cho, g, g);
    dist = _arb_vec_init(n * nb_z);
    t = _acb_vec_init(g);
    
    acb_theta_eld_cho(cho, tau, lp);
    for (k = 0; k < nb_z; k++)
    {
        acb_theta_dist_a0(dist + k * n, z + k * g, tau, lp);
    }

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&t[k]), state, prec);
        }
        res = acb_theta_ql_a0(r, t, z, dist, tau, guard, prec);
        guard += ACB_THETA_LOW_PREC;
    }

    /* Should actually be at 2tau; last duplication step, but solve the ql_step
       problems first */
   
    flint_randclear(state);
    arb_mat_clear(cho);
    _arb_vec_clear(dist, n * nb_z);
    _acb_vec_clear(t, g);
}
