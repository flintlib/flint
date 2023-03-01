/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

slong
acb_poly_valuation(const acb_poly_t poly)
{
    slong i, len = poly->length;

    for (i = 0; i < len; i++)
        if (!acb_is_zero(poly->coeffs + i))
            return i;

    return -1;
}

