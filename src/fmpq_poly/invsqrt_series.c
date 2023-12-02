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
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_invsqrt_series(fmpz * rpoly, fmpz_t rden,
                      const fmpz * poly, const fmpz_t den, slong len, slong n)
{
    slong m;
    fmpz * t, * u;
    fmpz_t tden, uden;

    if (n == 1)
    {
        fmpz_one(rpoly);
        fmpz_one(rden);
        return;
    }

    m = (n + 1) / 2;

    _fmpq_poly_invsqrt_series(rpoly, rden, poly, den, len, m);

    fmpz_init(tden);
    fmpz_init(uden);
    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);

    _fmpz_vec_zero(rpoly + m, n - m);

    _fmpq_poly_mul(t, tden, rpoly, rden, m, rpoly, rden, m);
    if (2*m - 1 < n)
        fmpz_zero(t + n - 1);

    _fmpq_poly_mullow(u, uden, t, tden, n, rpoly, rden, n, n);
    _fmpq_poly_mullow(t, tden, u, uden, n, poly, den, len, n);
    _fmpz_vec_neg(t + m, t + m, n - m);
    _fmpz_vec_zero(t, m);
    fmpz_mul_ui(tden, tden, UWORD(2));
    _fmpq_poly_canonicalise(t, tden, n);

    _fmpq_poly_add(rpoly, rden, rpoly, rden, m, t, tden, n);

    fmpz_clear(tden);
    fmpz_clear(uden);
    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
}

void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length < 1 || !fmpz_equal(poly->coeffs, poly->den))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_invsqrt_series). Constant term != 1.\n");
    }

    if (n < 1)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_invsqrt_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_invsqrt_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res); /* XXX: necessary? */
}
