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
acb_dirichlet_platt_lemma_B2(arb_t out, slong K, const arb_t h,
        const arb_t xi, slong prec)
{
    arb_t two, half, pi;
    arb_t x1, x2, x3, x4, x5;

    arb_init(two);
    arb_init(half);
    arb_init(pi);
    arb_init(x1);
    arb_init(x2);
    arb_init(x3);
    arb_init(x4);
    arb_init(x5);

    arb_set_ui(two, 2);
    arb_mul_2exp_si(half, two, -2);
    arb_const_pi(pi, prec);

    arb_set_si(x1, K + 5);
    arb_mul_2exp_si(x1, x1, -1);
    arb_pow(x1, two, x1, prec);

    arb_add_si(x2, half, K, prec);
    arb_pow(x2, pi, x2, prec);

    arb_pow_ui(x3, h, (ulong) (K + 1), prec);

    arb_pow_ui(x4, xi, (ulong) K, prec);

    arb_set_si(x5, K + 2);
    arb_mul_2exp_si(x5, x5, -1);
    arb_rgamma(x5, x5, prec);

    arb_mul(out, x1, x2, prec);
    arb_mul(out, out, x3, prec);
    arb_mul(out, out, x4, prec);
    arb_mul(out, out, x5, prec);

    arb_clear(two);
    arb_clear(half);
    arb_clear(pi);
    arb_clear(x1);
    arb_clear(x2);
    arb_clear(x3);
    arb_clear(x4);
    arb_clear(x5);
}
