/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq_poly.h"

void fmpq_poly_truncate(fmpq_poly_t poly, slong n)
{
    if (poly->length > n)
    {
        slong i;
        for (i = n; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);
        poly->length = n;
        fmpq_poly_canonicalise(poly);
    }
}
