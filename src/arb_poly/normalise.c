/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_normalise(arb_poly_t poly)
{
    slong i;

    for (i = poly->length - 1;
        (i >= 0) && arb_is_zero(poly->coeffs + i); i--);

    poly->length = i + 1;
}
