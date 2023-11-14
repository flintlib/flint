/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpq_poly.h"

static inline void __legendre_denom(fmpz_t den, ulong n)
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

void _fmpq_poly_legendre_p(fmpz * coeffs, fmpz_t den, ulong n)
{
    fmpz * r;
    int odd;
    slong k;
    ulong L;

    if (n == 0)
    {
        fmpz_one(coeffs);
        fmpz_one(den);
        return;
    }

    if (n == 1)
    {
        fmpz_zero(coeffs);
        fmpz_one(coeffs + 1);
        fmpz_one(den);
        return;
    }

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

void
fmpq_poly_legendre_p(fmpq_poly_t poly, ulong n)
{
    fmpq_poly_fit_length(poly, n + 1);
    _fmpq_poly_legendre_p(poly->coeffs, poly->den, n);
    _fmpq_poly_set_length(poly, n + 1);
}

