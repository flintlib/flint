/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

slong
arb_poly_valuation(const arb_poly_t poly)
{
    slong i, len = poly->length;

    for (i = 0; i < len; i++)
        if (!arb_is_zero(poly->coeffs + i))
            return i;

    return -1;
}

