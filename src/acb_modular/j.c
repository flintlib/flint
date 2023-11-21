/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_modular.h"

void
acb_modular_j(acb_t z, const acb_t tau, slong prec)
{
    psl2z_t g;
    arf_t one_minus_eps;
    acb_t tau_prime, t2, t3, t4, q;
    int real;

    if (!arb_is_positive(acb_imagref(tau)) || !arb_is_finite(acb_realref(tau)))
    {
        acb_indeterminate(z);
        return;
    }

    real = arb_is_int_2exp_si(acb_realref(tau), -1);

    psl2z_init(g);
    arf_init(one_minus_eps);
    acb_init(tau_prime);
    acb_init(t2);
    acb_init(t3);
    acb_init(t4);
    acb_init(q);

    arf_set_ui_2exp_si(one_minus_eps, 63, -6);
    acb_modular_fundamental_domain_approx(tau_prime, g, tau,
        one_minus_eps, prec);

    acb_exp_pi_i(q, tau_prime, prec);
    acb_modular_theta_const_sum(t2, t3, t4, q, prec);

    /* theta2 ^ 8 */
    acb_mul(t2, t2, t2, prec);
    acb_mul(t2, t2, t2, prec);
    acb_mul(t2, t2, q, prec);
    acb_mul(t2, t2, t2, prec);

    /* theta3 ^ 8 */
    acb_mul(t3, t3, t3, prec);
    acb_mul(t3, t3, t3, prec);
    acb_mul(t3, t3, t3, prec);

    /* theta4 ^ 8 */
    acb_mul(t4, t4, t4, prec);
    acb_mul(t4, t4, t4, prec);
    acb_mul(t4, t4, t4, prec);

    acb_mul(z, t2, t3, prec);
    acb_mul(z, z, t4, prec);

    acb_add(t2, t2, t3, prec);
    acb_add(t2, t2, t4, prec);
    acb_cube(t2, t2, prec);

    acb_div(z, t2, z, prec);
    acb_mul_2exp_si(z, z, 5);

    if (real)
        arb_zero(acb_imagref(z));

    psl2z_clear(g);
    arf_clear(one_minus_eps);
    acb_clear(tau_prime);
    acb_clear(t2);
    acb_clear(t3);
    acb_clear(t4);
    acb_clear(q);
}

