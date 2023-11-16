/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_evaluate_horner(ca_t y, ca_srcptr f, slong len,
                           const ca_t x, ca_ctx_t ctx)
{
    if (len == 0)
    {
        ca_zero(y, ctx);
    }
    else if (len == 1 || ca_check_is_zero(x, ctx) == T_TRUE)
    {
        ca_set(y, f, ctx);
    }
    else if (len == 2)
    {
        ca_mul(y, x, f + 1, ctx);
        ca_add(y, y, f + 0, ctx);
    }
    else
    {
        slong i = len - 1;
        ca_t t, u;

        ca_init(t, ctx);
        ca_init(u, ctx);
        ca_set(u, f + i, ctx);

        for (i = len - 2; i >= 0; i--)
        {
            ca_mul(t, u, x, ctx);
            ca_add(u, f + i, t, ctx);
        }

        ca_swap(y, u, ctx);

        ca_clear(t, ctx);
        ca_clear(u, ctx);
    }
}

void
ca_poly_evaluate_horner(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)
{
    _ca_poly_evaluate_horner(res, f->coeffs, f->length, a, ctx);
}
