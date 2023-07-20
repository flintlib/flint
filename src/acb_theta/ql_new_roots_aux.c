/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* In this function, guarantee nb_z >= 1 and z starts with 0 */

void acb_theta_ql_new_roots_aux(acb_ptr r, acb_ptr t, acb_srcptr z, slong nb_z,
    arb_srcptr dist, const acb_mat_t tau, slong nb_steps, slong guard, slong prec)
{    
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    flint_rand_t state;
    acb_ptr x;
    slong k, j;
    int res = 0;

    x = _acb_vec_init(2 * nb_z * g);
    flint_randinit(state);

    for (j = 0; (j < ACB_THETA_QL_TRY) && !res; j++)
    {
        /* Get z', 2z' picked at random in [0,2] */
        _acb_vec_zero(x, 2 * nb_z * g);
        for (k = 0; k < g; k++)
        {
            arb_urandom(acb_realref(&x[k]), state, prec);     
        }
        _acb_vec_scalar_mul_2exp_si(x, x, g, 1);
        _acb_vec_scalar_mul_2exp_si(x + g, x, g, 1);
        
        /* Set x */
        for (k = 1; k < nb_z; k++)
        {
            _acb_vec_add(x + k * 2 * g, x, z + k * g, g, prec);
            _acb_vec_add(x + k * 2 * g + g, x + g, z + k * g, g, prec);
        }

        /* Get roots */
        res = 1;
        for (k = 0; (k < 2 * nb_z) && res; k++)
        {
            res = acb_theta_ql_new_roots(r + k * nb_steps * n,
                x + k * g, dist + (k / 2) * n, tau, nb_steps, guard);
        }
    }
    
    _acb_vec_set(t, x, g);
    if (!res)
    {
        _acb_vec_indeterminate(r, 2 * nb_z * n * nb_steps);
        _acb_vec_indeterminate(t, g);
    }
    
    _acb_vec_clear(x, 2 * nb_z * g);
    flint_randclear(state);
}


