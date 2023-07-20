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
acb_theta_ql_new_roots(acb_ptr r, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t w;
    acb_ptr x;
    arb_t d;    
    slong hprec;
    slong k, a;
    int res = 1;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(g);
    arb_init(d);
   
    for (k = 0; k < nb_steps; k++)
    {
        acb_mat_scalar_mul_2exp_si(w, tau, k);
        _acb_vec_scalar_mul_2exp_si(x, z, g, k);

        for (a = 0; a < n; a++)
        {
            arb_const_log2(d, prec);
            arb_div(d, &dist[a], d, prec);
            arb_mul_2exp_si(d, d, k);
            hprec = prec + arf_get_si(arb_midref(d), ARF_RND_NEAR);
            acb_theta_naive_ind(&r[k * n + a], a << g, x, 1, w, hprec);

            flint_printf("(ql_new_roots) k = %wd, a = %wd, hprec = %wd, get:\n", k, a, hprec);
            acb_printd(&r[k * n + a], 10);
            flint_printf("\n");
            
            if (acb_contains_zero(&r[k * n + a]))
            {
                res = 0;
                break;
            }
        }
        if (res == 0)
        {
            break;
        }
    }

    acb_mat_clear(w);
    _acb_vec_clear(x, g);
    arb_clear(d);
    return res;
}
