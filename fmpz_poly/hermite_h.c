/*
    Copyright (C) 2011 Ralf Stephan
    Copyright (C) 2021 Mathieu Gouttenoire

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_hermite_h(fmpz * coeffs, ulong n)
{
    long k;
    
    if (n == 0)
    {
        fmpz_one(coeffs);
        return;
    }
    
    if (n == 1)
    {
        fmpz_zero(coeffs);
        fmpz_set_ui(coeffs + 1, 2);
        return;
    }
    
    for (k = n & 1; k < n; k += 2)
    {
        fmpz_zero(coeffs + k);
    }
    
    fmpz_one(coeffs + n);
    fmpz_mul_2exp(coeffs + n, coeffs + n, n);
    
    for (k = n - 2; k >= 0; k -= 2)
    {
        fmpz_mul2_uiui(coeffs + k, coeffs + k+2, k+1, k+2);
        fmpz_divexact_ui(coeffs + k, coeffs + k, (n - k) << 1);
        fmpz_neg(coeffs + k, coeffs + k);
    }
}

void
fmpz_poly_hermite_h(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_hermite_h(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n + 1);
}
