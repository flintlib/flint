/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "arb.h"
#include "bernoulli.h"

void
arb_zeta_ui_bernoulli(arb_t x, ulong n, slong prec)
{
    fmpq_t b;
    arb_t t, f;
    slong wp;

    if (n % 2)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    wp = prec + FLINT_BIT_COUNT(n) + 2;

    fmpq_init(b);
    arb_init(t);
    arb_init(f);

    bernoulli_fmpq_ui(b, n);
    arb_set_fmpq(x, b, wp);

    arb_const_pi(t, wp);
    arb_mul_2exp_si(t, t, 1);
    arb_pow_ui(t, t, n, wp);

    arb_fac_ui(f, n, wp);

    arb_div(t, t, f, wp);
    arb_mul(x, x, t, wp);
    arb_abs(x, x);
    arb_mul_2exp_si(x, x, -1);

    arb_clear(t);
    arb_clear(f);
    fmpq_clear(b);
}

