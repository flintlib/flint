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
acb_chebyshev_t2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)
{
    int i;

    acb_set_round(a, x, prec);
    acb_one(b);

    if (n <= 1)
    {
        if (n == 0)
            acb_swap(a, b);
        return;
    }

    for (i = FLINT_BIT_COUNT(n - 1) - 1; i >= 0; i--)
    {
        if (((n - 1) >> i) & 1)
        {
            acb_mul(b, b, a, prec);
            acb_mul_2exp_si(b, b, 1);
            acb_sub(b, b, x, prec);
            acb_mul(a, a, a, prec);
            acb_mul_2exp_si(a, a, 1);
            acb_sub_ui(a, a, 1, prec);
        }
        else
        {
            acb_mul(a, a, b, prec);
            acb_mul_2exp_si(a, a, 1);
            acb_sub(a, a, x, prec);
            acb_mul(b, b, b, prec);
            acb_mul_2exp_si(b, b, 1);
            acb_sub_ui(b, b, 1, prec);
        }
    }
}

