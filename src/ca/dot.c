/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca.h"

void
ca_dot(ca_t res, const ca_t initial, int subtract,
    ca_srcptr x, slong xstep, ca_srcptr y, slong ystep, slong len, ca_ctx_t ctx)
{
    slong i;
    ca_t t;

    if (len <= 0)
    {
        if (initial == NULL)
            ca_zero(res, ctx);
        else
            ca_set(res, initial, ctx);
        return;
    }

    ca_init(t, ctx);

    if (initial == NULL)
    {
        ca_mul(res, x, y, ctx);
    }
    else
    {
        if (subtract)
            ca_neg(res, initial, ctx);
        else
            ca_set(res, initial, ctx);

        ca_mul(t, x, y, ctx);
        ca_add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        ca_mul(t, x + i * xstep, y + i * ystep, ctx);
        ca_add(res, res, t, ctx);
    }

    if (subtract)
        ca_neg(res, res, ctx);

    ca_clear(t, ctx);
}
