/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

void
acb_dot_simple(acb_t res, const acb_t initial, int subtract,
    acb_srcptr x, slong xstep, acb_srcptr y, slong ystep, slong len, slong prec)
{
    slong i;

    if (len <= 0)
    {
        if (initial == NULL)
            acb_zero(res);
        else
            acb_set_round(res, initial, prec);
        return;
    }

    if (initial == NULL)
    {
        acb_mul(res, x, y, prec);
    }
    else
    {
        if (subtract)
            acb_neg(res, initial);
        else
            acb_set(res, initial);
        acb_addmul(res, x, y, prec);
    }

    for (i = 1; i < len; i++)
        acb_addmul(res, x + i * xstep, y + i * ystep, prec);

    if (subtract)
        acb_neg(res, res);
}
