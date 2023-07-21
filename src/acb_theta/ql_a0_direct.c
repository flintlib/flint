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
acb_theta_ql_a0_direct_init(acb_ptr r, acb_srcptr z, slong nb_z, arb_srcptr dist,
    const acb_mat_t tau, slong d, slong nb_steps, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    slong nb_a = 1 << (g - d);
    acb_mat_t w;
    acb_ptr x;
    arb_ptr new_dist;
    slong k, a;

    acb_mat_init(w, g, g);
    x = _acb_vec_init(nb_z * g);
    new_dist = _arb_vec_init(nb_z * n);

    acb_mat_scalar_mul_2exp_si(w, tau, nb_steps);
    _acb_vec_scalar_mul_2exp_si(x, z, nb_z * g, nb_steps);
    _arb_vec_scalar_mul_2exp_si(new_dist, dist, nb_z * n, nb_steps);

    for (k = 0; k < nb_z; k++)
    {
        for (a = 0; a < nb_a; a++)
        {
            /* Find out radius of ellipsoid to make */
        }
        if (d == 0)
        {
            /* Use naive algorithm */
        }
        else
        {
            /* Make ellipsoids */
        }
            
        
    }

    
    
}

/* In this function, assume nb_z >= 1 and z starts with 0 */
int acb_theta_ql_a0_direct(acb_ptr r, acb_srcptr z, slong nb_z,
    arb_srcptr dist, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t cho;
    acb_ptr roots;
    arb_ptr new_dist;
    slong* cuts;
    slong nb_cuts;
    slong d;
    slong nb_steps;
    slong k, j;
    int res = 1;

    arb_mat_init(cho, g, g);
    new_dist = _arb_vec_init(n);
    cuts = flint_malloc(g * sizeof(slong));

    acb_mat_get_imag(cho, tau);
    arb_mat_cho(cho, cho, prec);
    arb_mat_transpose(cho, cho);

    nb_cuts = acb_theta_ql_cuts(cuts, cho, prec);
    d = (nb_cuts > 0) ? cuts[nb_cuts - 1] : 0;
    nb_steps = acb_theta_ql_new_nb_steps(cho, d, prec);

    roots = _acb_vec_init(n * nb_z * nb_steps);

    for (k = 0; (k < nb_z) && res; k++)
    {
        res = acb_theta_ql_new_roots(roots + k * nb_steps * n, z + k * g,
            dist + k * n, tau, nb_steps, prec);
    }

    if (res)
    {
        res = acb_theta_ql_a0_direct_init(r, z, nb_z, dist, tau, d, nb_steps, prec);
    }

    if (res)
    {
        for (k = nb_steps - 1; k >= 0; k--)
        {
            for (j = 1; j < nb_z; j++)
            {
                _arb_vec_mul_2exp_si(new_dist, dist + j * n, n, k);
                acb_theta_ql_step(r + j * n, r + j * n, r,
                    roots + j * nb_steps * n + k * n, new_dist, g, prec);
            }
            _arb_vec_mul_2exp_si(new_dist, dist, n, k);
            acb_theta_ql_step(r, r, r, roots + k * n, new_dist, g, prec);
        }
    }

    arb_mat_clear(cho);
    _arb_vec_clear(new_dist, n);
    _acb_vec_clear(roots, n * nb_z * nb_steps);
    flint_free(cuts);
    return res;
}
