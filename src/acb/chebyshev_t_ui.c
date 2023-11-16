/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_chebyshev_t_ui(acb_t y, ulong n, const acb_t x, slong prec)
{
    int i, r;

    if (n <= 1)
    {
        if (n == 0)
            acb_one(y);
        else
            acb_set_round(y, x, prec);
        return;
    }

    r = flint_ctz(n);

    if ((n >> r) == 1)
    {
        acb_mul(y, x, x, prec);
        acb_mul_2exp_si(y, y, 1);
        acb_sub_ui(y, y, 1, prec);
        r -= 1;
    }
    else
    {
        /* we only need one value, so break out final iteration */
        acb_t t, u;
        acb_init(t);
        acb_init(u);
        acb_chebyshev_t2_ui(t, u, (n >> (r + 1)) + 1, x, prec);
        acb_mul(t, t, u, prec);
        acb_mul_2exp_si(t, t, 1);
        acb_sub(y, t, x, prec);
        acb_clear(t);
        acb_clear(u);
    }

    for (i = 0; i < r; i++)
    {
        acb_mul(y, y, y, prec);
        acb_mul_2exp_si(y, y, 1);
        acb_sub_ui(y, y, 1, prec);
    }
}

