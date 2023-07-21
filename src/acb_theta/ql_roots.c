/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static int
acb_theta_ql_roots_one(acb_ptr r, acb_srcptr z, arb_srcptr dist,
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
        arb_mul_2exp_si(d, &dist[a], k);

        for (a = 0; a < n; a++)
        {
            hprec = prec + acb_theta_ql_addprec(d);
            acb_theta_naive_ind(&r[k * n + a], a << g, x, 1, w, hprec);

            flint_printf("(ql_roots_one) k = %wd, a = %wd, hprec = %wd, get:\n", k, a, hprec);
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

int
acb_theta_ql_roots(acb_ptr r, acb_ptr t, acb_srcptr z, arb_srcptr dist,
    const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr x;
    slong k;
    int res = 1;

    if (_acb_vec_is_zero(t, g))
    {
        return acb_theta_ql_roots_one(r, z, dist, tau, nb_steps, guard, prec);
    }

    x = _acb_vec_init(g);
    
    for (k = 0; (k < 3) && res; k++)
    {
        _acb_vec_scalar_mul_si(x, t, g, k, prec);
        _acb_vec_add(x, x, z, g, prec);
        res = acb_theta_ql_roots_one(r + k * nb_steps * n, x, dist, tau,
            nb_steps, guard, prec);
    }

    _acb_vec_clear(x, g);
    return res;
}
