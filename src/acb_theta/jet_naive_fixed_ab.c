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
acb_theta_jet_naive_fixed_ab(acb_ptr dth, ulong ab, acb_srcptr z, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nb = acb_theta_jet_nb(ord, g);
    ulong a = ab >> g;
    ulong b = ab;
    acb_ptr v, w, new_z, aux;
    arb_ptr u;
    acb_t c, x;

    v = _acb_vec_init(g);
    w = _acb_vec_init(g);
    new_z = _acb_vec_init(g);
    aux = _acb_vec_init(nb);
    u = _arb_vec_init(g);
    acb_init(c);
    acb_init(x);

    acb_theta_char_get_acb(v, a, g);
    acb_theta_char_get_acb(w, b, g);
    acb_theta_char_get_arb(u, a, g);
    _arb_vec_scalar_mul_2exp_si(u, u, g, 1);

    /* Get jet at new_z */
    acb_mat_vector_mul_col(new_z, tau, v, prec);
    _acb_vec_add(new_z, new_z, w, g, prec);
    _acb_vec_add(new_z, new_z, z, g, prec);
    acb_theta_jet_naive_00(dth, new_z, tau, ord, prec);

    /* Get exponential factor */
    acb_mat_vector_mul_col(v, tau, v, prec);
    acb_theta_char_dot_acb(c, a, v, g, prec);
    _acb_vec_add(w, w, z, g, prec);
    acb_theta_char_dot_acb(x, a, w, g, prec);
    acb_mul_2exp_si(x, x, 1);
    acb_add(x, x, c, prec);
    acb_exp_pi_i(x, x, prec);

    /* Get other coefficients */
    acb_theta_jet_exp_pi_i(aux, u, ord, g, prec);
    _acb_vec_scalar_mul(aux, aux, nb, x, prec);
    acb_theta_jet_mul(dth, dth, aux, ord, g, prec);

    _acb_vec_clear(new_z, g);
    _acb_vec_clear(v, g);
    _acb_vec_clear(w, g);
    _acb_vec_clear(aux, nb);
    _arb_vec_clear(u, g);
    acb_clear(c);
    acb_clear(x);
}
