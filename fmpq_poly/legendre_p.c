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

    Copyright (C) 2016  Ralf Stephan

******************************************************************************/

#include "fmpq_poly.h"

void
_fmpq_poly_legendre_p(fmpz * coeffs, fmpz_t den, ulong n)
{
    fmpz_t c;
    ulong k = 0, kk = 0;

    if (n == 0)
    {
        fmpz_one(coeffs);
        return;
    }

    if (n == 1)
    {
        fmpz_one(coeffs + 1);
        return;
    }

    fmpz_init(c);
    fmpz_bin_uiui(c, 2*n, n);
    fmpz_zero(coeffs);
    fmpz_one(den);
    fmpz_set(coeffs + n, c);
    for (k = 1; k <= n/2; k++)
    {
        ++kk;
        fmpz_zero(coeffs + n - kk);

        ++kk;
        fmpz_mul2_uiui(c, c, n-2*k+1, n-2*k+2);
        fmpz_mul_ui(c, c, n-k+1);
        fmpz_divexact2_uiui(c, c, 2*n-2*k+1, 2*n-2*k+2);
        if (k>0)
            fmpz_divexact_ui(c, c, k);
        fmpz_neg(c, c);
        fmpz_set(coeffs + n - kk, c);
    }

    fmpz_one(c);
    fmpz_mul_2exp(c, c, n);
    _fmpq_poly_scalar_div_fmpz(coeffs, den, coeffs, den, n+1, c);

    fmpz_clear(c);
}

void
fmpq_poly_legendre_p(fmpq_poly_t poly, ulong n)
{
    fmpq_poly_fit_length(poly, n + 1);
    _fmpq_poly_legendre_p(poly->coeffs, poly->den, n);
    _fmpq_poly_set_length(poly, n + 1);
}
