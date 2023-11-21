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
acb_poly_fit_length(acb_poly_t poly, slong len)
{
    slong i;

    if (len > poly->alloc)
    {
        if (len < 2 * poly->alloc)
            len = 2 * poly->alloc;

        poly->coeffs = flint_realloc(poly->coeffs,
            len * sizeof(acb_struct));

        for (i = poly->alloc; i < len; i++)
            acb_init(poly->coeffs + i);

        poly->alloc = len;
    }
}
