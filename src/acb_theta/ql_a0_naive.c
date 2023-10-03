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
acb_theta_ql_a0_naive(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist0,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_t = !_acb_vec_is_zero(t, g);
    int has_z = !_acb_vec_is_zero(z, g);
    slong nb_t = (has_t ? 3 : 1);
    acb_ptr x, th;
    slong j, k;

    x = _acb_vec_init(g * nb_t);
    th = _acb_vec_init(nb_t);

    for (k = 0; k < nb_t; k++)
    {
        _acb_vec_scalar_mul_ui(x + k * g, t, g, k, prec);
    }
    for (k = 0; k < n; k++)
    {
        acb_theta_naive_fixed_ab(th, k << g, x, nb_t, tau,
            prec + acb_theta_dist_addprec(&dist0[k]));
        for (j = 0; j < nb_t; j++)
        {
            acb_set(&r[j * n + k], &th[j]);
        }
    }

    if (has_z)
    {
        for (k = 0; k < nb_t; k++)
        {
            _acb_vec_add(x + k * g, x + k * g, z, g, prec);
        }
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_fixed_ab(th, k << g, x, nb_t, tau,
                prec + acb_theta_dist_addprec(&dist[k]));
            for (j = 0; j < nb_t; j++)
            {
                acb_set(&r[nb_t * n + j * n + k], &th[j]);
            }
        }
    }

    _acb_vec_clear(x, g * nb_t);
    _acb_vec_clear(th, nb_t);
    return 1;
}
