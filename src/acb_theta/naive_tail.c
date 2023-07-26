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
acb_theta_naive_tail(arf_t bound, const arf_t R2, const arb_mat_t cho, slong ord, slong prec)
{
    arb_t res, temp;
    arb_t Rmod;
    slong g = arb_mat_nrows(cho);
    slong k;

    arb_init(res);
    arb_init(temp);
    arb_init(Rmod);

    /* Ensure assumptions R2\geq 4, R2\geq 2*ord are satisfied */
    arb_set_arf(Rmod, R2);
    arb_set_si(temp, FLINT_MAX(4, 2 * ord));
    arb_max(Rmod, Rmod, temp, prec);

    /* Evaluate 2^(2*g+2) R^(g-1 + 2*ord) exp(-R^2) \prod(1 + gamma_i^{-1}) */
    arb_one(res);
    arb_mul_2exp_si(res, res, 2 * g + 2);

    arb_sqrt(temp, Rmod, prec);
    arb_pow_ui(temp, temp, g - 1 + 2 * ord, prec);
    arb_mul(res, res, temp, prec);

    arb_neg(temp, Rmod);
    arb_exp(temp, temp, prec);
    arb_mul(res, res, temp, prec);

    for (k = 0; k < g; k++)
    {
        arb_inv(temp, arb_mat_entry(cho, k, k), prec);
        arb_add_si(temp, temp, 1, prec);
        arb_mul(res, res, temp, prec);
    }
    arb_get_ubound_arf(bound, res, prec);

    arb_clear(res);
    arb_clear(temp);
    arb_clear(Rmod);
}
