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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "arith.h"

void
_arith_chebyshev_t_polynomial(fmpz * coeffs, ulong n)
{
    long k, i, d, m;

    d = n % 2;

    fmpz_zero(coeffs);
    fmpz_set_ui(coeffs + d, d ? n : 1);
    if (n % 4 >= 2)
        fmpz_neg(coeffs + d, coeffs + d);

    m = n / 2;

    for (k = 1; k <= m; k++)
    {
        i = 2 * k + d;
        fmpz_mul2_uiui(coeffs + i, coeffs + i - 2, 4*(m-k+1), n+k-m-1);
        fmpz_divexact2_uiui(coeffs + i, coeffs + i, n+2*k-2*m-1, n+2*k-2*m);
        fmpz_neg(coeffs + i, coeffs + i);
        fmpz_zero(coeffs + i - 1);
    }
}

void
arith_chebyshev_t_polynomial(fmpz_poly_t poly, ulong n)
{
    if (n == 0)
    {
        fmpz_poly_set_ui(poly, 1UL);
        return;
    }

    fmpz_poly_fit_length(poly, n + 1);

    if (n == 1)
    {
        fmpz_zero(poly->coeffs);
        fmpz_one(poly->coeffs + 1);
    }
    else
        _arith_chebyshev_t_polynomial(poly->coeffs, n);

    _fmpz_poly_set_length(poly, n + 1);
}
