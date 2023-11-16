/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "acb.h"
#include "acb_modular.h"

void
acb_modular_eta(acb_t z, const acb_t tau, slong prec)
{
    psl2z_t g;
    fmpq_t t;
    arf_t one_minus_eps;
    acb_t tau_prime, eta, q, q24;

    psl2z_init(g);
    fmpq_init(t);
    arf_init(one_minus_eps);
    acb_init(tau_prime);
    acb_init(eta);
    acb_init(q);
    acb_init(q24);

    arf_set_ui_2exp_si(one_minus_eps, 63, -6);
    acb_modular_fundamental_domain_approx(tau_prime, g, tau,
        one_minus_eps, prec);

    acb_div_ui(q24, tau_prime, 12, prec);
    acb_exp_pi_i(q24, q24, prec);
    acb_pow_ui(q, q24, 24, prec);

    acb_modular_eta_sum(eta, q, prec);
    acb_mul(eta, eta, q24, prec);

    /* epsilon^-1 */
    fmpq_set_si(t, -acb_modular_epsilon_arg(g), 12);
    arb_sin_cos_pi_fmpq(acb_imagref(q), acb_realref(q), t, prec);

    acb_mul(eta, eta, q, prec);

    /* (c*tau+d)^(-1/2) */
    if (!fmpz_is_zero(&g->c))
    {
        acb_mul_fmpz(q, tau, &g->c, prec);
        acb_add_fmpz(q, q, &g->d, prec);
        acb_rsqrt(q, q, prec);
        acb_mul(eta, eta, q, prec);
    }

    acb_set(z, eta);

    psl2z_clear(g);
    fmpq_clear(t);
    arf_clear(one_minus_eps);
    acb_clear(tau_prime);
    acb_clear(eta);
    acb_clear(q);
    acb_clear(q24);
}

