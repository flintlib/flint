/*
    Copyright (C) 2021 Albin Ahlbäck, Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

/* Up to this order the coefficients fit inside small fmpz */
#define SMALL_BOUND (FLINT_BITS == 64 ? 20 : 12)

static void
_fmpz_vec_powers(fmpz * res, ulong n, slong len)
{
    slong ix;

    if (len != 0)
        fmpz_set_ui(res + 0, n == 0);

    for (ix = 1; ix < len; ix += 2)
    {
        fmpz_set_ui(res + ix, ix);
        fmpz_pow_ui(res + ix, res + ix, n);
    }

    for (ix = 2; ix < len; ix += 2)
        fmpz_mul_2exp(res + ix, res + ix / 2, n);
}

/* Only valid for when 1 <= len <= n / 2 */
static void
_fmpz_vec_binomials(fmpz * res, ulong n, slong len)
{
    slong ix, jx;

    fmpz_one(res + 0);
    fmpz_set_ui(res + 1, n);

    for (ix = 2, jx = 1; ix < len; ix++, jx++)
    {
        fmpz_mul_ui(res + ix, res + jx, n - jx);
        fmpz_divexact_si(res + ix, res + ix, ix);
    }
}

void
_fmpz_poly_eulerian_polynomial_series(fmpz * res, ulong n)
{
    slong ix, m;
    fmpz * tmp;

    m = (n + 1) / 2;
    tmp = _fmpz_vec_init(2 * m + 1);
    _fmpz_vec_binomials(tmp, n + 1, m);
    for (ix = 1; ix < m; ix += 2)
        fmpz_neg(tmp + ix, tmp + ix);
    _fmpz_vec_powers(tmp + m, n, m + 1);
    _fmpz_poly_mullow(res, tmp, m, tmp + m + 1, m, m);
    _fmpz_vec_clear(tmp, 2 * m + 1);
}

void
_fmpz_poly_eulerian_polynomial_rec(fmpz * res, ulong n)
{
    slong ix, jx;

    fmpz_one(res);
    for (ix = 1; ix <= FLINT_MIN(n / 2, SMALL_BOUND / 2); ix++)
        _fmpz_demote(res + ix);

    for (ix = 3; ix <= n; ix++)
    {
        if (ix <= SMALL_BOUND)
        {
            if (ix & 1)
                res[ix / 2] = (ix + 1) * res[ix / 2 - 1];
            for (jx = ix / 2 - 1; jx >= 1; jx--)
                res[jx] = (ix - jx) * res[jx - 1] + (jx + 1) * res[jx];
        }
        else
        {
            if (ix & 1)
                fmpz_mul_ui(res + ix / 2, res + ix / 2 - 1, ix + 1);
            for (jx = ix / 2 - 1; jx >= 1; jx--)
            {
                fmpz_mul_ui(res + jx, res + jx, jx + 1);
                fmpz_addmul_ui(res + jx, res + (jx - 1), ix - jx);
            }
        }
    }
}

void
_fmpz_poly_eulerian_polynomial(fmpz * res, ulong n)
{
    slong ix;

    if (n < 32)
        _fmpz_poly_eulerian_polynomial_rec(res, n);
    else
        _fmpz_poly_eulerian_polynomial_series(res, n);

    for (ix = 0; ix < n / 2; ix++)
        fmpz_set(res + (n - ix - 1), res + ix);
}

void
fmpz_poly_eulerian_polynomial(fmpz_poly_t poly, ulong n)
{
    /* A_0, A_1 = 1 */
    if (n <= 1)
    {
        fmpz_poly_one(poly);
        return;
    }

    fmpz_poly_fit_length(poly, n);
    _fmpz_poly_eulerian_polynomial(poly->coeffs, n);
    _fmpz_poly_set_length(poly, n);
}
