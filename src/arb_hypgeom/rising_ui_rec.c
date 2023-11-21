/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

void
arb_hypgeom_rising_ui_rec(arb_t res, const arb_t x, ulong n, slong prec)
{
    if (n <= 1)
    {
        if (n == 0)
            arb_one(res);
        else
            arb_set_round(res, x, prec);
        return;
    }

    if (n == 2 && prec <= 1024)
    {
        if (res != x)
            arb_set(res, x);
        arb_addmul(res, x, x, prec);
        return;
    }

    if ((prec < 512 && n <= 20) || (n <= FLINT_MIN(80, 6000 / prec)))
    {
        arb_hypgeom_rising_ui_forward(res, x, n, prec);
    }
    else
    {
        if (n >= 20 && arb_bits(x) < prec / 8)
            arb_hypgeom_rising_ui_bs(res, x, n, prec);
        else
            arb_hypgeom_rising_ui_rs(res, x, n, 0, prec);
    }
}

