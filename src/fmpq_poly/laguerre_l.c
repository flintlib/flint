/*
    Copyright (C) 2016  Ralf Stephan

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"

void _fmpq_poly_laguerre_l(fmpz * coeffs, fmpz_t den, ulong n)
{
    fmpz_t c;
    ulong k;

    if (n == 0)
    {
        fmpz_one(coeffs);
        fmpz_one(den);
        return;
    }

    if (n == 1)
    {
        fmpz_one(coeffs);
        fmpz_one(coeffs + 1);
        fmpz_neg(coeffs + 1, coeffs + 1);
        fmpz_one(den);
        return;
    }

    fmpz_init(c);
    fmpz_one(c);
    if (n%2 == 1)
        fmpz_neg(c, c);
    fmpz_set(coeffs + n, c);

    for (k = 0; k < n; k++)
    {
        fmpz_mul2_uiui(c, c, n-k, n-k);
        fmpz_divexact_ui(c, c, k+1);
        fmpz_neg(c, c);
        fmpz_set(coeffs + n - k - 1, c);
    }
    fmpz_set(den, coeffs);
    fmpz_clear(c);
}

void
fmpq_poly_laguerre_l(fmpq_poly_t poly, ulong n)
{
    fmpq_poly_fit_length(poly, n + 1);
    _fmpq_poly_laguerre_l(poly->coeffs, poly->den, n);
    _fmpq_poly_set_length(poly, n + 1);
}

