/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

void gr_poly_init(gr_poly_t poly, gr_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
gr_poly_init2(gr_poly_t poly, slong len, gr_ctx_t ctx)
{
    gr_poly_init(poly, ctx);
    gr_poly_fit_length(poly, len, ctx);
}
