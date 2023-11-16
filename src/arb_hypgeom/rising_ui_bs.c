/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

static void
bsplit(arb_t y, const arb_t x, ulong a, ulong b, slong prec)
{
    if (b - a <= 16)
    {
        if (a == 0)
        {
            arb_hypgeom_rising_ui_forward(y, x, b, prec);
        }
        else
        {
            arb_add_ui(y, x, a, prec);
            arb_hypgeom_rising_ui_forward(y, y, b - a, prec);
        }
    }
    else
    {
        arb_t t, u;
        ulong m = a + (b - a) / 2;

        arb_init(t);
        arb_init(u);

        bsplit(t, x, a, m, prec);
        bsplit(u, x, m, b, prec);

        arb_mul(y, t, u, prec);

        arb_clear(t);
        arb_clear(u);
    }
}

void
arb_hypgeom_rising_ui_bs(arb_t res, const arb_t x, ulong n, slong prec)
{
    if (n <= 1)
    {
        if (n == 0)
            arb_one(res);
        else
            arb_set_round(res, x, prec);
        return;
    }

    {
        arb_t t;
        slong wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        arb_init(t);
        bsplit(t, x, 0, n, wp);
        arb_set_round(res, t, prec);
        arb_clear(t);
    }
}

