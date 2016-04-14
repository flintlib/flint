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

void _fmpq_poly_gegenbauer_c(fmpz * coeffs, fmpz_t den, ulong n, fmpq_t a)
{
    fmpq_t c;
    fmpz_t t, p, nu, de;
    ulong k;
    slong kk;

    if (n == 0)
    {
        fmpz_one(coeffs);
        fmpz_one(den);
        return;
    }

    if (n == 1)
    {
        fmpq_init(c);
        fmpz_zero(coeffs);
        fmpq_mul_2exp(c, a, 1);
        fmpq_canonicalise(c);
        fmpz_set(coeffs + 1, fmpq_numref(c));
        fmpz_set(den, fmpq_denref(c));
        fmpq_clear(c);
        return;
    }

    fmpz_init(t);
    fmpz_init(p);
    fmpz_init(nu);
    fmpz_init(de);

    fmpz_set(nu, fmpq_numref(a));
    fmpz_set(de, fmpq_denref(a));
    fmpz_pow_ui(den, de, n);
    fmpz_fac_ui(t, n);
    fmpz_mul(den, den, t);
    fmpz_fac_ui(p, n/2);
    fmpz_divexact(p, t, p);

    if (n%2)
        fmpz_mul_2exp(p, p, 1);
    if (n&2)
        fmpz_neg(p, p);

    for (k = 0; k < n-n/2; k++)
    {
        fmpz_mul(p, p, nu);
        fmpz_add(nu, nu, de);
    }

    fmpz_pow_ui(t, de, n/2);
    fmpz_mul(p, p, t);
    fmpz_set(coeffs + (n%2), p);

    for (kk = n/2 - 1; kk >= 0; --kk)
    {
        fmpz_mul(p, p, nu);
        fmpz_divexact(p, p, de);
        fmpz_divexact2_uiui(p, p, n-2*kk-1, n-2*kk);
        fmpz_mul_ui(p, p, kk+1);
        fmpz_mul_2exp(p, p, 2);
        fmpz_neg(p, p);
        fmpz_set(coeffs + n - 2*kk, p);
        fmpz_add(nu, nu, de);
    }

    fmpz_clear(t);
    fmpz_clear(p);
    fmpz_clear(nu);
    fmpz_clear(de);
}

void
fmpq_poly_gegenbauer_c(fmpq_poly_t poly, ulong n, fmpq_t a)
{
    fmpq_poly_fit_length(poly, n + 1);
    _fmpq_poly_gegenbauer_c(poly->coeffs, poly->den, n, a);
    _fmpq_poly_set_length(poly, n + 1);
}

