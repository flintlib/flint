/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, slong len1,
                              const fmpz * poly2, slong len2)
{
    _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, 0, len1 + len2 - 1);
}

void
fmpz_poly_mul_KS(fmpz_poly_t res,
                 const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const slong len1 = poly1->length;
    const slong len2 = poly2->length;
    const slong rlen = len1 + len2 - 1;

    if (len1 == 0 || len2 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    fmpz_poly_fit_length(res, rlen);
    _fmpz_poly_mul_KS(res->coeffs, poly1->coeffs, len1, poly2->coeffs, len2);
    _fmpz_poly_set_length(res, rlen);
}

