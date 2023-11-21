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
acb_poly_init(acb_poly_t poly)
{
    poly->coeffs = NULL;
    poly->length = 0;
    poly->alloc = 0;
}

void
acb_poly_init2(acb_poly_t poly, slong len)
{
    acb_poly_init(poly);
    acb_poly_fit_length(poly, len);
}
