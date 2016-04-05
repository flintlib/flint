/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Ralf Stephan

******************************************************************************/

#include "fmpz_poly.h"

void
_fmpz_poly_hermite_h(fmpz * coeffs, ulong n)
{
    fmpz_t c, two;
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
    fmpz_init(two);
    fmpz_set_ui(two, 2);
    fmpz_pow_ui(c, two, n);
    fmpz_clear(two);

    while (1)
    {
        fmpz_set(coeffs + n, c);
        fmpz_divexact_si(c, c, -4);
        fmpz_mul_ui(c, c, n);
        if (--n == 0)
            break;

        fmpz_zero(coeffs + n);
        fmpz_mul_ui(c, c, n);
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
