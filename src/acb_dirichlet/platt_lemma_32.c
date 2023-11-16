/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* out = 2*pi^(5/4)*exp(1/(8*h^2) - (t0^2)/(2*h^2) - pi*x) */
void
acb_dirichlet_platt_lemma_32(arb_t out, const arb_t h, const arb_t t0,
        const arb_t x, slong prec)
{
    arb_t pi, one_fourth;
    arb_t x1, x2;

    arb_init(pi);
    arb_init(one_fourth);
    arb_init(x1);
    arb_init(x2);

    arb_const_pi(pi, prec);
    arb_set_d(one_fourth, 0.25);

    arb_set_d(x1, 1.25);
    arb_pow(x1, pi, x1, prec);
    arb_mul_2exp_si(x1, x1, 1);

    arb_sqr(x2, t0, prec);
    arb_sub(x2, x2, one_fourth, prec);
    arb_div(x2, x2, h, prec);
    arb_div(x2, x2, h, prec);
    arb_mul_2exp_si(x2, x2, -1);

    arb_mul(out, pi, x, prec);
    arb_add(out, out, x2, prec);
    arb_neg(out, out);
    arb_exp(out, out, prec);
    arb_mul(out, out, x1, prec);

    arb_clear(pi);
    arb_clear(one_fourth);
    arb_clear(x1);
    arb_clear(x2);
}
