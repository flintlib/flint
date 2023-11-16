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
ca_factor_get_ca(ca_t res, const ca_factor_t fac, ca_ctx_t ctx)
{
    slong i, len;
    ca_t t;

    len = fac->length;

    if (len == 0)
    {
        ca_one(res, ctx);
    }
    else if (len == 1)
    {
        ca_pow(res, fac->base, fac->exp, ctx);
    }
    else
    {
        ca_init(t, ctx);
        ca_pow(res, fac->base, fac->exp, ctx);

        for (i = 1; i < len; i++)
        {
            ca_pow(t, fac->base + i, fac->exp + i, ctx);
            ca_mul(res, res, t, ctx);
        }

        ca_clear(t, ctx);
    }
}

