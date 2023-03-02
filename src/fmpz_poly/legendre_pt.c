/*
    Copyright (C) 2016  Ralf Stephan

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

void
_fmpz_poly_legendre_pt(fmpz * coeffs, ulong n)
{
    fmpz_t c;
    ulong k, facnp, facnm;

    if (n == 0)
    {
        fmpz_one(coeffs);
        return;
    }

    if (n == 1)
    {
        fmpz_set_si(coeffs, -1);
        fmpz_set_ui(coeffs + 1, 2);
        return;
    }

    fmpz_init(c);
    fmpz_one(c);
    if (n%2 == 1)
        fmpz_neg(c, c);
    facnp = n;
    facnm = n+1;
    fmpz_set(coeffs, c);

    for (k = 1; k<=n; k++)
    {
        ++facnp;
        --facnm;
        fmpz_mul2_uiui(c, c, facnp, facnm);
        fmpz_divexact2_uiui(c, c, k, k);
        fmpz_neg(c, c);
        fmpz_set(coeffs + k, c);
    }

    fmpz_clear(c);
}

void
fmpz_poly_legendre_pt(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_legendre_pt(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n + 1);
}
