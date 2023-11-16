/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_platt_lemma_B1(arb_t out,
        slong sigma, const arb_t t0, const arb_t h, const fmpz_t J, slong prec)
{
    arb_t pi, C, x1, x2, x3, Ja;

    if (sigma % 2 == 0 || sigma < 3)
    {
        arb_zero_pm_inf(out);
        return;
    }

    arb_init(pi);
    arb_init(C);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(Ja);

    arb_const_pi(pi, prec);
    acb_dirichlet_platt_c_bound(C, sigma, t0, h, 0, prec);
    arb_set_fmpz(Ja, J);

    arb_set_si(x1, 2*sigma - 1);
    arb_div(x1, x1, h, prec);
    arb_sqr(x1, x1, prec);
    arb_mul_2exp_si(x1, x1, -3);
    arb_exp(x1, x1, prec);

    arb_set_si(x2, 1 - 2*sigma);
    arb_mul_2exp_si(x2, x2, -2);
    arb_pow(x2, pi, x2, prec);

    arb_set_si(x3, 1 - sigma);
    arb_pow(x3, Ja, x3, prec);
    arb_div_si(x3, x3, sigma - 1, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x3, prec);
    arb_mul(out, out, C, prec);

    arb_clear(pi);
    arb_clear(C);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(Ja);
}
