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
acb_theta_jet_bounds_2(arb_t c, arb_t rho, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t mat;
    acb_mat_t w;
    arb_ptr y;
    arb_t t;
    arf_t rad;
    slong j, k;

    arb_mat_init(mat, g, g);
    acb_mat_init(w, g, g);
    y = _arb_vec_init(g);
    arb_init(t);
    arf_init(rad);

    acb_mat_get_imag(mat, tau);

    /* Get lower bound on radius around tau */
    arb_mat_spd_radius(rad, mat, prec);
    arf_mul_2exp_si(rad, rad, -1);
    arb_set_arf(rho, rad);

    /* Set w to matrix with larger error */
    acb_mat_set(w, tau);
    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
        {
            acb_add_error_arf(acb_mat_entry(w, j, k), rad);
        }
    }

    /* Get upper bound on exponential sum */
    acb_theta_eld_cho(mat, w, prec);
    arb_one(c);
    for (j = 0; j < g; j++)
    {
        arb_sqrt(t, arb_mat_entry(mat, j, j), prec);
        arb_inv(t, t, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_add_si(t, t, 1, prec);
        arb_mul(c, c, t, prec);
    }
    arb_mul_2exp_si(c, c, g);

    /* Multiply by exponential factor */
    acb_mat_get_imag(mat, w);
    arb_mat_inv(mat, mat, prec);
    arb_const_pi(t, prec);
    arb_mat_scalar_mul_arb(mat, mat, t, prec);
    _acb_vec_get_imag(y, z, g);
    for (j = 0; j < g; j++)
    {
        arb_add_error(&y[j], rho);
    }
    arb_mat_bilinear_form(t, mat, y, y, prec);
    arb_exp(t, t, prec);
    arb_mul(c, c, t, prec);

    arb_mat_clear(mat);
    acb_mat_clear(w);
    _arb_vec_clear(y, g);
    arb_clear(t);
    arf_clear(rad);
}
