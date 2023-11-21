/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"

void
acb_hypgeom_rising_ui_rec(acb_t res, const acb_t x, ulong n, slong prec)
{
    if (n <= 1)
    {
        if (n == 0)
            acb_one(res);
        else
            acb_set_round(res, x, prec);
        return;
    }

    if (arb_is_zero(acb_imagref(x)))
    {
        arb_hypgeom_rising_ui_rec(acb_realref(res), acb_realref(x), n, prec);
        arb_zero(acb_imagref(res));
        return;
    }

    if (n == 2 && prec <= 1024)
    {
        if (res != x)
            acb_set(res, x);
        acb_addmul(res, x, x, prec);
        return;
    }

    if (n <= 5 && prec <= 512)
    {
        acb_hypgeom_rising_ui_forward(res, x, n, prec);
    }
    else
    {
        if (n >= 20 && acb_bits(x) < prec / 8)
            acb_hypgeom_rising_ui_bs(res, x, n, prec);
        else
            acb_hypgeom_rising_ui_rs(res, x, n, 0, prec);
    }
}

