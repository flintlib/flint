/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

void
_gr_poly_set_length(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    if (poly->length > len)
    {
        GR_MUST_SUCCEED(_gr_vec_zero(GR_ENTRY(poly->coeffs, len, ctx->sizeof_elem), poly->length - len, ctx));
    }

    poly->length = len;
}
