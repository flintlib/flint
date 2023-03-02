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
arb_chebyshev_t_ui(arb_t y, ulong n, const arb_t x, slong prec)
{
    int i, r;

    if (n <= 1)
    {
        if (n == 0)
            arb_one(y);
        else
            arb_set_round(y, x, prec);
        return;
    }

    count_trailing_zeros(r, n);

    if ((n >> r) == 1)
    {
        arb_mul(y, x, x, prec);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
        r -= 1;
    }
    else
    {
        /* we only need one value, so break out final iteration */
        arb_t t, u;
        arb_init(t);
        arb_init(u);
        arb_chebyshev_t2_ui(t, u, (n >> (r + 1)) + 1, x, prec);
        arb_mul(t, t, u, prec);
        arb_mul_2exp_si(t, t, 1);
        arb_sub(y, t, x, prec);
        arb_clear(t);
        arb_clear(u);
    }

    for (i = 0; i < r; i++)
    {
        arb_mul(y, y, y, prec);
        arb_mul_2exp_si(y, y, 1);
        arb_sub_ui(y, y, 1, prec);
    }
}

