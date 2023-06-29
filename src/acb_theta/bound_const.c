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
acb_theta_bound_const(arf_t rad, arf_t bound, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    arb_mat_t im;
    acb_mat_t pert;
    arb_t lambda;
    slong j, k;

    arb_mat_init(im, g, g);
    acb_mat_init(pert, g, g);
    arb_init(lambda);

    acb_mat_get_imag(im, tau);

    /* Get lower bound on radius around tau */
    arb_mat_spd_neighborhood(rad, im, prec);
    arf_mul_2exp_si(rad, rad, -1);

    /* Get upper bound for exponential sum */
    acb_mat_set(pert, tau);
    for (j = 0; j < g; j++)
    {
        for (k = 0; k < g; k++)
            acb_add_error_arf(acb_mat_entry(pert, j, k), rad);
    }
    acb_mat_get_imag(im, pert);
    arb_mat_spd_eig_lbound_arf(arb_midref(lambda), im, prec);
    arb_sqrt(lambda, lambda, prec);
    arb_inv(lambda, lambda, prec);
    arb_add_si(lambda, lambda, 1, prec);
    arb_pow_ui(lambda, lambda, g, prec);
    arb_get_ubound_arf(bound, lambda, prec);

    arb_mat_clear(im);
    acb_mat_clear(pert);
    arb_clear(lambda);
}
