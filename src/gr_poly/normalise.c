/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

void
_gr_poly_normalise(gr_poly_t poly, gr_ctx_t ctx)
{
    slong i, sz;
    truth_t eq;

    i = poly->length - 1;
    sz = ctx->sizeof_elem;

    while (i >= 0)
    {
        eq = gr_is_zero(GR_ENTRY(poly->coeffs, i, sz), ctx);

        if (eq == T_TRUE)
        {
            GR_MUST_SUCCEED(gr_zero(GR_ENTRY(poly->coeffs, i, sz), ctx));
            i--;
        }
        else
        {
            break;
        }
    }

    poly->length = i + 1;
}
