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
acb_theta_naive_fixed_a(acb_ptr th, ulong a, acb_srcptr z, slong nb_z,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr new_z;
    acb_ptr v, w;
    acb_t c, x;
    slong k, b;

    new_z = _acb_vec_init(nb_z * g);
    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    acb_init(c);
    acb_init(x);

    acb_theta_char_get_acb(v, a, g);
    acb_mat_vector_mul_col(v, tau, v, prec); /* tau.a/2 */
    for (k = 0; k < nb_z; k++)
    {
        _acb_vec_add(new_z + k * g, z + k * g, v, g, prec);
    }

    acb_theta_naive_0b(th, new_z, nb_z, tau, prec);

    acb_theta_char_dot_acb(c, a, v, g, prec);
    for (k = 0; k < nb_z; k++)
    {
        for (b = 0; b < n; b++)
        {
            acb_theta_char_get_acb(w, b, g);
            _acb_vec_add(w, w, z + k * g, g, prec);
            acb_theta_char_dot_acb(x, a, w, g, prec);
            acb_mul_2exp_si(x, x, 1);
            acb_add(x, x, c, prec);
            acb_exp_pi_i(x, x, prec);
            acb_mul(&th[k * n + b], &th[k * n + b], x, prec);
        }
    }

    _acb_vec_clear(new_z, nb_z * g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    acb_clear(c);
    acb_clear(x);
}
