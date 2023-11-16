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
gr_poly_fit_length(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    slong alloc = poly->alloc;

    if (len > alloc)
    {
        slong sz = ctx->sizeof_elem;

        if (len < 2 * alloc)
            len = 2 * alloc;

        poly->coeffs = flint_realloc(poly->coeffs, len * sz);
        _gr_vec_init(GR_ENTRY(poly->coeffs, poly->alloc, sz), len - alloc, ctx);
        poly->alloc = len;
    }
}
