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

void _arb_increment_fast(arb_t x, slong prec);

void
acb_hypgeom_rising_ui_forward(acb_t res, const acb_t x, ulong n, slong prec)
{
    acb_t t;
    ulong k;
    slong wp;

    if (n <= 1)
    {
        if (n == 0)
            acb_one(res);
        else
            acb_set_round(res, x, prec);
        return;
    }

    wp = prec + FLINT_BIT_COUNT(n);

    acb_init(t);

    acb_add_ui(t, x, 1, wp);
    acb_mul(res, x, t, (n == 2) ? prec : wp);

    for (k = 2; k < n; k++)
    {
        _arb_increment_fast(acb_realref(t), wp);
        acb_mul(res, res, t, k == (n - 1) ? prec : wp);
    }

    acb_clear(t);
}

