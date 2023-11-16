/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_elliptic.h"
#include "acb_modular.h"

/* eta1 = (-1/6) theta1'''(0) / theta1'(0) */
static void
eta1(acb_t res, acb_t theta1, const acb_t tau, slong prec)
{
    acb_ptr theta;
    acb_t z;
    acb_init(z);
    theta = _acb_vec_init(16);
    acb_modular_theta_jet(theta,
        theta + 4, theta + 8, theta + 12, z, tau, 4, prec);
    if (theta1 != NULL)
        acb_set(theta1, theta + 1);
    acb_div(res, theta + 3, theta + 1, prec);
    acb_neg(res, res);
    _acb_vec_clear(theta, 16);
    acb_clear(z);
}

/* zeta(z) = 2 eta1 z + theta1'(z,q) / theta1(z,q) */
void
acb_elliptic_zeta(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_ptr t;
    int real;

    real = acb_is_real(z) && arb_is_int_2exp_si(acb_realref(tau), -1) &&
                arb_is_positive(acb_imagref(tau));

    t = _acb_vec_init(8);

    acb_modular_theta_jet(t, t + 2, t + 4, t + 6, z, tau, 2, prec);
    eta1(t + 2, NULL, tau, prec);
    acb_mul_2exp_si(t + 2, t + 2, 1);
    acb_mul(t + 2, t + 2, z, prec);
    acb_div(t, t + 1, t, prec);
    acb_add(res, t, t + 2, prec);
    if (real)
        arb_zero(acb_imagref(res));

    _acb_vec_clear(t, 8);
}

/* sigma(z) = exp(eta1 z^2) theta1(z) / theta1'(0) */
void
acb_elliptic_sigma(acb_t res, const acb_t z, const acb_t tau, slong prec)
{
    acb_ptr t;
    int real;

    real = acb_is_real(z) && arb_is_int_2exp_si(acb_realref(tau), -1) &&
                arb_is_positive(acb_imagref(tau));

    t = _acb_vec_init(8);

    acb_modular_theta_jet(t, t + 2, t + 4, t + 6, z, tau, 2, prec);
    eta1(t + 2, t + 3, tau, prec);
    acb_mul(t + 4, z, z, prec);
    acb_mul(t + 2, t + 2, t + 4, prec);
    acb_exp(t + 2, t + 2, prec);
    acb_div(t, t, t + 3, prec);
    acb_mul(res, t, t + 2, prec);
    if (real)
        arb_zero(acb_imagref(res));

    _acb_vec_clear(t, 8);
}

