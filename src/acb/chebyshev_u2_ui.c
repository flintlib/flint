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
acb_chebyshev_u2_ui(acb_t a, acb_t b, ulong n, const acb_t x, slong prec)
{
    int i;
    acb_t t, u;

    if (n == 0)
    {
        acb_one(a);
        acb_zero(b);
        return;
    }

    acb_set_round(a, x, prec);
    acb_mul_2exp_si(a, a, 1);
    acb_one(b);

    if (n == 1)
        return;

    acb_init(t);
    acb_init(u);

    for (i = FLINT_BIT_COUNT(n) - 2; i >= 0; i--)
    {
        acb_add(t, a, b, prec);
        acb_sub(u, a, b, prec);

        if ((n >> i) & 1)
        {
            acb_submul(b, x, a, prec);
            acb_mul(a, a, b, prec);
            acb_neg(a, a);
            acb_mul_2exp_si(a, a, 1);
            acb_mul(b, t, u, prec);
        }
        else
        {
            acb_submul(a, x, b, prec);
            acb_mul(b, a, b, prec);
            acb_mul_2exp_si(b, b, 1);
            acb_mul(a, t, u, prec);
        }
    }

    acb_clear(t);
    acb_clear(u);
}

