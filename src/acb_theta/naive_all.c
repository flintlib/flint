/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_naive_all(acb_ptr th, acb_srcptr z, slong nb_z, const acb_mat_t tau,
                    slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr all_z, ata, v;
    acb_t c;
    slong a, b, d, k;
    
    all_z = _acb_vec_init(g * n * nb_z);
    ata = _acb_vec_init(n);
    v = _acb_vec_init(g);
    acb_init(c);
    
    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_acb(v, a, g);
        acb_mat_vector_mul_col(v, tau, v, prec);
        for (k = 0; k < nb_z; k++)
        {
            _acb_vec_add(all_z + k * g * n + a * g, z + k * g, v, g, prec);
        }
        acb_char_dot_acb(&ata[a], a, v, g, prec);
    }

    acb_theta_naive_0b(th, all_z, n * nb_z, tau, prec);

    for (a = 0; a < n; a++)
    {
        /* Factors depending on z, not on b */
        for (k = 0; k < nb_z; k++)
        {
            acb_theta_char_dot_acb(c, a, z + k * g, g, prec);
            acb_mul_2exp_si(c, c, 1);
            acb_add(c, c, &ata[a], prec);
            acb_exp_pi_i(c, c, prec);
            _acb_vec_scalar_mul_acb(th + k * n * n + a * n,
                th + k * n * n + a * n, n, c, prec);
        }
        /* Factors depending on b, not on z */
        for (b = 0; b < n; k++)
        {
            d = acb_theta_char_dot(a, b, g);
            for (k = 0; k < nb_z; k++)
            {
                acb_mul_powi(&th[k * n * n + a * n + b],
                    &th[k * n * n + a * n + b], d);
            }
        }
    }

    _acb_vec_clear(all_z, g * n * nb_z);
    _acb_vec_clear(ata, n);
    _acb_vec_clear(v, g);
    acb_clear(c);
}
