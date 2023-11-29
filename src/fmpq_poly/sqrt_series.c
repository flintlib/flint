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
#include "fmpz_poly.h"
#include "fmpq_poly.h"


void
_fmpq_poly_sqrt_series(fmpz * rpoly, fmpz_t rden,
                      const fmpz * poly, const fmpz_t den, slong len, slong n)
{
    fmpz * t;
    fmpz_t tden;
    t = _fmpz_vec_init(n);
    fmpz_init(tden);
    _fmpq_poly_invsqrt_series(t, tden, poly, den, len, n);
    _fmpq_poly_mullow(rpoly, rden, t, tden, n, poly, den, len, n);
    _fmpz_vec_clear(t, n);
    fmpz_clear(tden);
}

void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length < 1 || !fmpz_equal(poly->coeffs, poly->den))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sqrt_series). Constant term != 1.\n");
    }

    if (n < 1)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_sqrt_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_sqrt_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    fmpq_poly_canonicalise(res);
}
