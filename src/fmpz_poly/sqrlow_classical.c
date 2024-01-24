/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    Assumes len > 0 and 0 < n <= 2 * len - 1.
 */
void _fmpz_poly_sqrlow_classical(fmpz *res, const fmpz *op, slong len, slong n)
{
    slong i, start, stop;

    len = FLINT_MIN(len, n);

    if (n == 1)
    {
        fmpz_mul(res, op, op);
        return;
    }

    fmpz_mul(res, op, op);
    fmpz_mul(res + 1, op, op + 1);
    fmpz_mul_2exp(res + 1, res + 1, 1);

    for (i = 2; i < FLINT_MIN(n, 2 * len - 3); i++)
    {
        start = FLINT_MAX(0, i - len + 1);
        stop = FLINT_MIN(len - 1, (i + 1) / 2 - 1);

        _fmpz_vec_dot_general(res + i, NULL, 0, op + start, op + i - stop, 1, stop - start + 1);
        fmpz_mul_2exp(res + i, res + i, 1);

        if (i % 2 == 0 && i / 2 < len)
            fmpz_addmul(res + i, op + i / 2, op + i / 2);
    }

    if (len > 2 && n >= 2 * len - 2)
    {
        fmpz_mul(res + 2 * len - 3, op + len - 1, op + len - 2);
        fmpz_mul_2exp(res + 2 * len - 3, res + 2 * len - 3, 1);
    }

    if (n >= 2 * len - 1)
        fmpz_mul(res + 2 * len - 2, op + len - 1, op + len - 1);
}

void
fmpz_poly_sqrlow_classical(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
{
    slong len = poly->length;

    if (len == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    n = FLINT_MIN(2 * len - 1, n);

    if (res == poly)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_sqrlow_classical(t->coeffs, poly->coeffs, len, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_sqrlow_classical(res->coeffs, poly->coeffs, len, n);
    }

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
