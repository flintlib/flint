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
acb_theta_ql_a0(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int has_z = !_acb_vec_is_zero(z, g);
    int has_t = !_acb_vec_is_zero(t, g);
    slong nb_z = (has_z ? 2 : 1);
    slong nb_t = (has_t ? 3 : 1);
    acb_ptr x, th;
    arb_ptr dist0;
    int res;

    x = _acb_vec_init(g);
    th = _acb_vec_init(nb_z * nb_t * n);
    dist0 = _arb_vec_init(n);

    if (has_z)
    {
        acb_theta_dist_a0(dist0, x, tau, ACB_THETA_LOW_PREC);
    }
    else
    {
        _arb_vec_set(dist0, dist, n);
    }
    res = acb_theta_ql_a0_steps(th, t, z, dist, dist0, tau, guard, prec,
        &acb_theta_ql_a0);

    if (has_z)
    {
        _acb_vec_set(r, th + nb_t * n, n * nb_t);
    }
    else
    {
        _acb_vec_set(r, th, n * nb_t);
    }

    _acb_vec_clear(x, g);
    _acb_vec_clear(th, nb_z * nb_t * n);
    _arb_vec_clear(dist0, n);
    return res;
}
