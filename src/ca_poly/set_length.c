/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_poly.h"

void
_ca_poly_set_length(ca_poly_t poly, slong len, ca_ctx_t ctx)
{
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            ca_zero(poly->coeffs + i, ctx);
    }

    poly->length = len;
}
