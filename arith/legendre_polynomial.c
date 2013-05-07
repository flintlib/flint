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
#include "fmpq_poly.h"
#include "arith.h"

static __inline__ void __legendre_denom(fmpz_t den, ulong n)
{
    ulong d, k;
    d = k = n >> 1;

    while (k)
    {
        k >>= 1;
        d += k;
    }

    fmpz_one(den);
    fmpz_mul_2exp(den, den, d);
}

void _arith_legendre_polynomial(fmpz * coeffs, fmpz_t den, ulong n)
{
    fmpz * r;
    int odd;
    long k;
    ulong L;

    L = n / 2;
    odd = n % 2;

    r = coeffs + odd;

    __legendre_denom(den, n);

    fmpz_bin_uiui(r, n, L);
    fmpz_mul(r, r, den);
    if (odd)
        fmpz_mul_ui(r, r, L + 1);
    fmpz_fdiv_q_2exp(r, r, 2*L);
    if (L % 2)
        fmpz_neg(r, r);

    for (k = 1; k <= L; k++)
    {
        fmpz_mul2_uiui(r + 2, r, L + 1 - k, 2*k + 2*L - 1 + 2*odd);
        fmpz_divexact2_uiui(r + 2, r + 2, k, 2*k - 1 + 2*odd);
        fmpz_neg(r + 2, r + 2);
        r += 2;
    }

    for (k = 1 - odd; k < n; k += 2)
        fmpz_zero(coeffs + k);
}

void arith_legendre_polynomial(fmpq_poly_t poly, ulong n)
{
    if (n == 0)
    {
        fmpq_poly_set_ui(poly, 1UL);
        return;
    }

    fmpq_poly_fit_length(poly, n + 1);

    if (n == 1)
    {
        fmpz_zero(poly->coeffs);
        fmpz_one(poly->coeffs + 1);
        fmpz_one(poly->den);
    }
    else
        _arith_legendre_polynomial(poly->coeffs, poly->den, n);

    _fmpq_poly_set_length(poly, n + 1);
}
