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
acb_modular_lambda(acb_t r, const acb_t tau, slong prec)
{
    psl2z_t g;
    arf_t one_minus_eps;
    acb_t tau_prime, q;
    acb_struct thetas[4];
    int R[4], S[4], C;
    int Rsum, qpower;

    psl2z_init(g);
    arf_init(one_minus_eps);
    acb_init(tau_prime);
    acb_init(q);
    acb_init(thetas + 0);
    acb_init(thetas + 1);
    acb_init(thetas + 2);
    acb_init(thetas + 3);

    arf_set_ui_2exp_si(one_minus_eps, 63, -6);
    acb_modular_fundamental_domain_approx(tau_prime, g, tau,
        one_minus_eps, prec);

    acb_modular_theta_transform(R, S, &C, g);

    acb_exp_pi_i(q, tau_prime, prec);
    acb_modular_theta_const_sum(thetas + 1, thetas + 2, thetas + 3, q, prec);
    acb_zero(thetas + 0);

    /* divide the transformation factors */
    Rsum = 4 * (R[1] - R[2]);
    /* possible factor [q^(+/- 1/4)]^4 needed for theta_1^4 or theta_2^4 */
    qpower = (S[1] == 0 || S[1] == 1) - (S[2] == 0 || S[2] == 1);

    acb_div(r, thetas + S[1], thetas + S[2], prec);
    acb_mul(r, r, r, prec);
    acb_mul(r, r, r, prec);

    if ((Rsum & 7) == 4)
        acb_neg(r, r);

    if (qpower == 1)
        acb_mul(r, r, q, prec);
    else if (qpower == -1)
        acb_div(r, r, q, prec);

    psl2z_clear(g);
    arf_clear(one_minus_eps);
    acb_clear(tau_prime);
    acb_clear(q);
    acb_clear(thetas + 0);
    acb_clear(thetas + 1);
    acb_clear(thetas + 2);
    acb_clear(thetas + 3);
}

