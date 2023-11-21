/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_hypgeom.h"

static void
bsplit(acb_t y, const acb_t x, ulong a, ulong b, slong prec)
{
    if (b - a <= 4)
    {
        if (a == 0)
        {
            acb_hypgeom_rising_ui_forward(y, x, b, prec);
        }
        else
        {
            acb_add_ui(y, x, a, prec);
            acb_hypgeom_rising_ui_forward(y, y, b - a, prec);
        }
    }
    else
    {
        acb_t t, u;
        ulong m = a + (b - a) / 2;

        acb_init(t);
        acb_init(u);

        bsplit(t, x, a, m, prec);
        bsplit(u, x, m, b, prec);

        acb_mul(y, t, u, prec);

        acb_clear(t);
        acb_clear(u);
    }
}

void
acb_hypgeom_rising_ui_bs(acb_t res, const acb_t x, ulong n, slong prec)
{
    if (n <= 1)
    {
        if (n == 0)
            acb_one(res);
        else
            acb_set_round(res, x, prec);
        return;
    }

    {
        acb_t t;
        slong wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        acb_init(t);
        bsplit(t, x, 0, n, wp);
        acb_set_round(res, t, prec);
        acb_clear(t);
    }
}

