/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_chebyshev_u(fmpz * coeffs, ulong n)
{
    slong k, i, d, m;

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

    d = n % 2;

    fmpz_zero(coeffs);
    fmpz_set_ui(coeffs + d, d ? n + 1 : 1);
    if (n % 4 >= 2)
        fmpz_neg(coeffs + d, coeffs + d);

    m = n / 2;

    for (k = 1; k <= m; k++)
    {
        i = 2 * k + d;
        fmpz_mul2_uiui(coeffs + i, coeffs + i - 2, 4*(m-k+1), n+k-m);
        fmpz_divexact2_uiui(coeffs + i, coeffs + i, n+2*k-2*m-1, n+2*k-2*m);
        fmpz_neg(coeffs + i, coeffs + i);
        fmpz_zero(coeffs + i - 1);
    }
}

void
fmpz_poly_chebyshev_u(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_chebyshev_u(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n + 1);
}

