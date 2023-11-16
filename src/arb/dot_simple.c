/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
arb_dot_simple(arb_t res, const arb_t initial, int subtract,
    arb_srcptr x, slong xstep, arb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            arb_zero(res);
        else
            arb_set_round(res, initial, prec);
        return;
    }

    if (initial == NULL)
    {
        arb_mul(res, x, y, prec);
    }
    else
    {
        if (subtract)
            arb_neg(res, initial);
        else
            arb_set(res, initial);
        arb_addmul(res, x, y, prec);
    }

    for (i = 1; i < len; i++)
        arb_addmul(res, x + i * xstep, y + i * ystep, prec);

    if (subtract)
        arb_neg(res, res);
}
