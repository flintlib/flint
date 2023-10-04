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
acb_theta_ql_a0_naive(acb_ptr th, acb_srcptr t, acb_srcptr z, arb_srcptr d0,
    arb_srcptr d, const acb_mat_t tau, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    int hast = !_acb_vec_is_zero(t, g);
    int hasz = !_acb_vec_is_zero(z, g);
    slong nbt = (hast ? 3 : 1);
    acb_ptr x, aux;
    slong j, k;

    x = _acb_vec_init(g * nbt);
    aux = _acb_vec_init(nbt);

    for (k = 0; k < nbt; k++)
    {
        _acb_vec_scalar_mul_ui(x + k * g, t, g, k, prec);
    }
    for (k = 0; k < n; k++)
    {
        acb_theta_naive_fixed_ab(aux, k << g, x, nbt, tau,
            prec + acb_theta_dist_addprec(&d0[k]));
        for (j = 0; j < nbt; j++)
        {
            acb_set(&th[j * n + k], &aux[j]);
        }
    }

    if (hasz)
    {
        for (k = 0; k < nbt; k++)
        {
            _acb_vec_add(x + k * g, x + k * g, z, g, prec);
        }
        for (k = 0; k < n; k++)
        {
            acb_theta_naive_fixed_ab(aux, k << g, x, nbt, tau,
                prec + acb_theta_dist_addprec(&d[k]));
            for (j = 0; j < nbt; j++)
            {
                acb_set(&th[nbt * n + j * n + k], &aux[j]);
            }
        }
    }

    _acb_vec_clear(x, g * nbt);
    _acb_vec_clear(aux, nbt);
    return 1;
}
