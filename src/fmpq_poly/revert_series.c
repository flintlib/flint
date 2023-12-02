/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_revert_series(fmpz * Qinv, fmpz_t den,
        const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    if (fmpz_is_one(Qden) && (n > 1) && fmpz_is_pm1(Q + 1))
    {
        _fmpz_poly_revert_series(Qinv, Q, Qlen, n);
        fmpz_one(den);
        return;
    }

    _fmpq_poly_revert_series_lagrange_fast(Qinv, den, Q, Qden, Qlen, n);
}

void
fmpq_poly_revert_series(fmpq_poly_t res,
            const fmpq_poly_t poly, slong n)
{
    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        flint_throw(FLINT_ERROR, "(fmpq_poly_revert_series): "
                "Input must have zero constant term and nonzero coefficient of x^1.\n");
    }

    if (n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_revert_series(res->coeffs,
                res->den, poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series(t->coeffs,
                t->den, poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
