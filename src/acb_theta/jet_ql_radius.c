/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_theta.h"

void
acb_theta_jet_ql_radius(arf_t eps, arf_t err, const arb_t c, const arb_t rho,
    slong ord, slong g, slong prec)
{
    slong lp = ACB_THETA_LOW_PREC;
    slong b = ord + 1;
    arb_t x, y;

    arb_init(x);
    arb_init(y);

    /* Set x to min of (1/2g)^(1/b)*rho, (2^(-prec)/2cg)^(1/b)*rho^(2b-1)/b */
    arb_set_si(x, 2 * g);
    arb_inv(x, x, lp);
    arb_root_ui(x, x, b, lp);
    arb_mul(x, x, rho, lp);

    arb_pow_ui(y, rho, 2 * b - 1, prec);
    arb_mul_2exp_si(y, y, -prec);
    arb_div(y, y, c, lp);
    arb_div_si(y, y, 2 * g, lp);
    arb_root_ui(y, y, b, lp);

    arb_min(x, x, y, lp);
    arb_get_lbound_arf(eps, x, lp);

    arf_one(err);
    arf_mul_2exp_si(err, err, -prec);

    arb_clear(x);
    arb_clear(y);
}
