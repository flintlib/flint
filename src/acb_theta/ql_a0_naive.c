/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int
acb_theta_ql_a0_naive(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    slong nb_z = (has_t ? 3 : 1);
    acb_ptr x, th;
    slong j, k;

    x = _acb_vec_init(nb_z * g);
    th = _acb_vec_init(nb_z);

    _acb_vec_set(x, z, g);
    if (has_t)
    {
        _acb_vec_set(x + g, t, g);
        _acb_vec_add(x + g, x + g, z, g, prec);
        _acb_vec_scalar_mul_2exp_si(x + 2 * g, t, g, 1);
        _acb_vec_add(x + 2 * g, x + 2 * g, z, g, prec);
    }

    for (k = 0; k < n; k++)
    {
        acb_theta_naive_ind(th, k << g, x, nb_z, tau,
            prec + acb_theta_dist_addprec(&dist[k]));
        for (j = 0; j < nb_z; j++)
        {
            acb_set(&r[j * n + k], &th[j]);
        }
    }

    _acb_vec_clear(x, nb_z * g);
    _acb_vec_clear(th, nb_z);
    return 1;
}
