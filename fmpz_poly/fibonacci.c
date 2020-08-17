/*
    Copyright (C) 2016 Shivin Srivastava

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void _fmpz_poly_fibonacci(fmpz * coeffs, ulong n)
{
    fmpz * r;
    int even;
    slong k;
    ulong L;

    if (n == 0) return;

    if (n == 1)
    {
        fmpz_one(coeffs);
        return;
    }

    L = n / 2;
    even = 1 - (n % 2);

    /* set the first two coefficients of poly depending parity of n */
    if (even)
    {
        fmpz_zero(coeffs);
        fmpz_one(coeffs + 1);
        fmpz_mul_ui(coeffs + 1, coeffs + 1, L);
    }
    else 
    {
        fmpz_one(coeffs);
        fmpz_zero(coeffs + 1);
    }

    fmpz_one(coeffs + n - 1);
    
    r = coeffs + even;
    r += 2;

    /* calculate the coefficients of the polynomial*/
    for (k = 2 + even; k < n - 2; k += 2)
    {
        fmpz_mul2_uiui(r, r - 2, L + k / 2, L + k / 2 - k + 1);
        fmpz_divexact2_uiui(r, r, k, k - 1);
        r += 2;
    }

    /* set the alternate coefficients to 0 again depending on the parity*/
    for (k = 1 + even; k < n; k += 2)
    {
        fmpz_zero(coeffs + k);
    }
}

void fmpz_poly_fibonacci(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n);
    _fmpz_poly_fibonacci(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n);
}

