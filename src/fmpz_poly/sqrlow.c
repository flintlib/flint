/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson

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


void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, slong len, slong n)
{
    if (len == 1)
    {
        fmpz_mul(res, poly, poly);
        return;
    }

    _fmpz_poly_mulmid(res, poly, len, poly, len, 0, n);
}

void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    const slong len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        fmpz_poly_sqrlow(t, poly, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
        return;
    }

    n = FLINT_MIN(2 * len - 1, n);

    fmpz_poly_fit_length(res, n);
    _fmpz_poly_sqrlow(res->coeffs, poly->coeffs, len, n);
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
