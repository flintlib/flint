/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "padic_poly.h"

void padic_poly_truncate(padic_poly_t poly, slong n, const fmpz_t p)
{
    if (poly->length > n)
    {
        slong i;

        for (i = n; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = n;
        _padic_poly_normalise(poly);
        padic_poly_canonicalise(poly, p);
    }
}
