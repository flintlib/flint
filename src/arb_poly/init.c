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
arb_poly_init(arb_poly_t poly)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
arb_poly_init2(arb_poly_t poly, slong len)
{
    arb_poly_init(poly);
    arb_poly_fit_length(poly, len);
}
