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
arb_chebyshev_t2_ui(arb_t a, arb_t b, ulong n, const arb_t x, slong prec)
{
    int i;

    arb_set_round(a, x, prec);
    arb_one(b);

    if (n <= 1)
    {
        if (n == 0)
            arb_swap(a, b);
        return;
    }

    for (i = FLINT_BIT_COUNT(n - 1) - 1; i >= 0; i--)
    {
        if (((n - 1) >> i) & 1)
        {
            arb_mul(b, b, a, prec);
            arb_mul_2exp_si(b, b, 1);
            arb_sub(b, b, x, prec);
            arb_mul(a, a, a, prec);
            arb_mul_2exp_si(a, a, 1);
            arb_sub_ui(a, a, 1, prec);
        }
        else
        {
            arb_mul(a, a, b, prec);
            arb_mul_2exp_si(a, a, 1);
            arb_sub(a, a, x, prec);
            arb_mul(b, b, b, prec);
            arb_mul_2exp_si(b, b, 1);
            arb_sub_ui(b, b, 1, prec);
        }
    }
}

