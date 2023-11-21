/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_set_length(acb_poly_t poly, slong len)
{
    slong i;

    if (poly->length > len)
    {
        for (i = len; i < poly->length; i++)
            acb_zero(poly->coeffs + i);
    }

    poly->length = len;
}
