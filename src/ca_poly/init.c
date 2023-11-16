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
ca_poly_init(ca_poly_t poly, ca_ctx_t ctx)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
ca_poly_init2(ca_poly_t poly, slong len, ca_ctx_t ctx)
{
    ca_poly_init(poly, ctx);
    ca_poly_fit_length(poly, len, ctx);
}
