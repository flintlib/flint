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
worker_d(acb_ptr r, acb_srcptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr x, th;
    int res;

    x = _acb_vec_init(2 * g);
    th = _acb_vec_init(6 * n);

    _acb_vec_set(x + g, z, g);
    acb_theta_ql_a0_with_t(th, t, x, 2,
        

    _acb_vec_clear(x, 2 * g);
    _acb_vec_clear(th, 6 * n);
}

int
acb_theta_ql_a0_with_t(acb_ptr r, acb_srcptr t, acb_srcptr z, slong nb_z,
    arb_srcptr dist, const acb_mat_t tau, slong guard, slong prec)
{
    
}
