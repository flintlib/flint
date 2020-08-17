/*
    Copyright (C) 2011 Ralf Stephan

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
    fmpz_t c;
    ulong fac = 1;

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

    fmpz_init(c);
    fmpz_one(c);
    fmpz_mul_2exp(c, c, n);

    while (1)
    {
        fmpz_set(coeffs + n, c);
        if (--n == 0)
            break;

        fmpz_zero(coeffs + n);
        fmpz_neg(c, c);
        fmpz_mul2_uiui(c, c, n+1, n);
        fmpz_fdiv_q_2exp(c, c, 2);
        fmpz_divexact_ui(c, c, fac);
        ++fac;
        if (--n == 0)
        {
            fmpz_set(coeffs, c);
            break;
        }
    }

    fmpz_clear(c);
}

void
fmpz_poly_hermite_h(fmpz_poly_t poly, ulong n)
{
    fmpz_poly_fit_length(poly, n + 1);
    _fmpz_poly_hermite_h(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n + 1);
}
