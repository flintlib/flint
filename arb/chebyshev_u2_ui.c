/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_chebyshev_u2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)
{
    int i;
    arb_t t, u;

    if (n == 0)
    {
        arb_one(a);
        arb_zero(b);
        return;
    }

    arb_set_round(a, x, prec);
    arb_mul_2exp_si(a, a, 1);
    arb_one(b);

    if (n == 1)
        return;

    arb_init(t);
    arb_init(u);

    for (i = FLINT_BIT_COUNT(n) - 2; i >= 0; i--)
    {
        arb_add(t, a, b, prec);
        arb_sub(u, a, b, prec);

        if ((n >> i) & 1)
        {
            arb_submul(b, x, a, prec);
            arb_mul(a, a, b, prec);
            arb_neg(a, a);
            arb_mul_2exp_si(a, a, 1);
            arb_mul(b, t, u, prec);
        }
        else
        {
            arb_submul(a, x, b, prec);
            arb_mul(b, a, b, prec);
            arb_mul_2exp_si(b, b, 1);
            arb_mul(a, t, u, prec);
        }
    }

    arb_clear(t);
    arb_clear(u);
}

