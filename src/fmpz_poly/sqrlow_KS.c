/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_sqrlow_KS(fmpz * res, const fmpz * poly, slong len, slong n)
{
    _fmpz_poly_mullow_KS(res, poly, len, poly, len, n);
}

void fmpz_poly_sqrlow_KS(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    const slong len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    n = FLINT_MIN(n, 2 * len - 1);
    fmpz_poly_fit_length(res, n);
    _fmpz_poly_sqrlow_KS(res->coeffs, poly->coeffs, len, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
