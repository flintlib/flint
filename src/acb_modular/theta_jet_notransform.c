/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

void
acb_modular_theta_jet_notransform(acb_ptr theta1, acb_ptr theta2,
    acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau,
    slong len, slong prec)
{
    acb_t q, q4, w;
    int w_is_unit;

    acb_init(q);
    acb_init(q4);
    acb_init(w);

    /* compute q_{1/4}, q */
    acb_mul_2exp_si(q4, tau, -2);
    acb_exp_pi_i(q4, q4, prec);
    acb_pow_ui(q, q4, 4, prec);

    /* compute w */
    acb_exp_pi_i(w, z, prec);
    w_is_unit = arb_is_zero(acb_imagref(z));

    /* evaluate theta functions */
    acb_modular_theta_sum(theta1, theta2, theta3, theta4,
        w, w_is_unit, q, len, prec);
    _acb_vec_scalar_mul(theta1, theta1, len, q4, prec);
    _acb_vec_scalar_mul(theta2, theta2, len, q4, prec);

    acb_clear(q);
    acb_clear(q4);
    acb_clear(w);
}

