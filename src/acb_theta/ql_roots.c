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
_acb_vec_entry_contains_zero(acb_srcptr v, slong len)
{
    slong k;
    
    for (k = 0; k < len; k++)
    {
        if (acb_contains_zero(&v[k]))
        {
            return 1;
        }
    }
    return 0;
}

slong
acb_theta_ql_roots(acb_ptr r, acb_srcptr z, slong nb_z, const acb_mat_t tau,
    slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong gap = acb_theta_ql_roots_max_gap(g);
    acb_mat_t w;
    acb_ptr x;
    slong hprec;
    slong k, j;
    int fail;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g * nb_z);
    
    hprec = 10;
    fail = 0;
    for (k = 0; (k < nb_steps) && !fail; k++)
    {
        acb_mat_scalar_mul_2exp_si(w, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, z, g * nb_z, k);
        fail = 1;
        for (j = 0; j < gap; j++)
        {
            acb_theta_naive_a0(r + k * n * nb_z, z, nb_z, w, hprec);
            if (!_acb_vec_entry_contains_zero(r + k * n * nb_z, n * nb_z))
            {
                fail = 0;
                break;
            }
            hprec *= 2;
        }
    }
    hprec += prec + 4 * g * n_clog(nb_steps, 2);

    acb_mat_clear(w);
    _acb_vec_clear(x, g * nb_z);
    return (fail ? -1 : hprec);
}
