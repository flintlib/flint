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
ca_set_d_d(ca_t res, double x, double y, ca_ctx_t ctx)
{
    if (y == 0.0)
    {
        ca_set_d(res, x, ctx);
    }
    else
    {
        ca_t t;
        ca_init(t, ctx);
        ca_set_d(t, y, ctx);
        ca_i(res, ctx);
        ca_mul(res, res, t, ctx);
        ca_set_d(t, x, ctx);
        ca_add(res, res, t, ctx);
        ca_clear(t, ctx);
    }
}
