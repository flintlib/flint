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
acb_theta_naive_tail(arb_t res, const arf_t R2, const arb_mat_t C, slong ord)
{
    slong g = arb_mat_nrows(C);
    slong lp = ACB_THETA_LOW_PREC;
    arb_t t, Rm;
    slong k;

    arb_init(t);
    arb_init(Rm);

    /* Ensure assumptions R2\geq 4, R2\geq 2*ord are satisfied */
    arb_set_arf(Rm, R2);
    arb_set_si(t, FLINT_MAX(4, 2 * ord));
    arb_max(Rm, Rm, t, lp);

    /* Evaluate 2^(2*g+2) R^(g - 1 + ord) exp(-R^2) \prod(1 + gamma_i^{-1}) */
    arb_one(res);
    arb_mul_2exp_si(res, res, 2 * g + 2);

    arb_sqrt(t, Rm, lp);
    arb_pow_ui(t, t, g - 1 + ord, lp);
    arb_mul(res, res, t, lp);

    arb_neg(t, Rm);
    arb_exp(t, t, lp);
    arb_mul(res, res, t, lp);

    for (k = 0; k < g; k++)
    {
        arb_inv(t, arb_mat_entry(C, k, k), lp);
        arb_add_si(t, t, 1, lp);
        arb_mul(res, res, t, lp);
    }

    arb_clear(t);
    arb_clear(Rm);
}
