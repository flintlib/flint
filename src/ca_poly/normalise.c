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
_ca_poly_normalise(ca_poly_t poly, ca_ctx_t ctx)
{
    slong i;

    i = poly->length - 1;

    while (i >= 0 && (ca_check_is_zero(poly->coeffs + i, ctx) == T_TRUE))
    {
        ca_zero(poly->coeffs + i, ctx);
        i--;
    }

    poly->length = i + 1;
}
