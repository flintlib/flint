/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_naive_fixed_ab(acb_ptr th, ulong ab, acb_srcptr zs, slong nb,
    const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    ulong a = ab >> g;
    ulong b = ab;
    acb_ptr new_zs;
    acb_ptr v, w;
    acb_t c, x;
    slong k;

    new_zs = _acb_vec_init(nb * g);
    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    acb_init(c);
    acb_init(x);

    acb_theta_char_get_acb(v, a, g);
    acb_mat_vector_mul_col(v, tau, v, prec); /* tau.a/2 */
    acb_theta_char_get_acb(w, b, g);
    _acb_vec_add(w, v, w, g, prec);
    for (k = 0; k < nb; k++)
    {
        _acb_vec_add(new_zs + k * g, zs + k * g, w, g, prec);
    }

    acb_theta_naive_00(th, new_zs, nb, tau, prec);

    acb_theta_char_dot_acb(c, a, v, g, prec);
    for (k = 0; k < nb; k++)
    {
        acb_theta_char_get_acb(w, b, g);
        _acb_vec_add(w, w, zs + k * g, g, prec);
        acb_theta_char_dot_acb(x, a, w, g, prec);
        acb_mul_2exp_si(x, x, 1);
        acb_add(x, x, c, prec);
        acb_exp_pi_i(x, x, prec);
        acb_mul(&th[k], &th[k], x, prec);
    }

    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    acb_clear(c);
    acb_clear(x);
}
