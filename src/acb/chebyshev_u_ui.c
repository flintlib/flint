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
acb_chebyshev_u_ui(acb_t y, ulong n, const acb_t x, slong prec)
{
    acb_t a, b;

    if (n <= 1)
    {
        if (n == 0)
        {
            acb_one(y);
        }
        else
        {
            acb_set_round(y, x, prec);
            acb_mul_2exp_si(y, y, 1);
        }
        return;
    }

    acb_init(a);
    acb_init(b);

    acb_chebyshev_u2_ui(a, b, n / 2, x, prec);

    if (n % 2 == 0)
    {
        acb_add(y, a, b, prec);
        acb_sub(b, a, b, prec);
        acb_mul(y, y, b, prec);
    }
    else
    {
        acb_submul(b, a, x, prec);
        acb_mul(y, a, b, prec);
        acb_mul_2exp_si(y, y, 1);
        acb_neg(y, y);
    }

    acb_clear(a);
    acb_clear(b);
}

