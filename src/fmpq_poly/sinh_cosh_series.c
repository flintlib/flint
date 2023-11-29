/*
    Copyright (C) 2016, 2021 Fredrik Johansson

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
_fmpq_poly_sinh_cosh_series(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden,
    slong Alen, slong n)
{
    fmpz * t;
    fmpz_t tden;

    t = _fmpz_vec_init(n);
    fmpz_init(tden);

    /* sinh(x) = (exp(x)-exp(-x))/2 */
    /* cosh(x) = (exp(x)+exp(-x))/2 = sinh(x) + exp(-x) */
    _fmpq_poly_exp_expinv_series(S, Sden, t, tden, A, Aden, Alen, n);
    _fmpq_poly_sub(S, Sden, S, Sden, n, t, tden, n);
    _fmpq_poly_scalar_div_ui(S, Sden, S, Sden, n, 2);
    _fmpq_poly_add(C, Cden, S, Sden, n, t, tden, n);

    _fmpz_vec_clear(t, n);
    fmpz_clear(tden);
}

void
fmpq_poly_sinh_cosh_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_zero(res2);
        return;
    }

    if (fmpq_poly_is_zero(poly) || n == 1)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_one(res2);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sinh_cosh_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res1, n);
    fmpq_poly_fit_length(res2, n);
    _fmpq_poly_sinh_cosh_series(res1->coeffs, res1->den,
        res2->coeffs, res2->den, poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res1, n);
    _fmpq_poly_normalise(res1);
    _fmpq_poly_set_length(res2, n);
    _fmpq_poly_normalise(res2);
}
