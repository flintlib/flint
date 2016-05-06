/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

/* Requires 2*min(Qlen,n) + n - 1 < 3n coefficients of scratch space in W */
static void
_fmpz_poly_inv_series_basecase_rev(fmpz * Qinv, fmpz * W,
    const fmpz * Q, slong Qlen, slong n)
{
    slong Wlen;
    fmpz *Qrev;

    Qlen = FLINT_MIN(Qlen, n);
    Wlen = n + Qlen - 1;
    Qrev = W + Wlen;

    _fmpz_poly_reverse(Qrev, Q, Qlen, Qlen);
    _fmpz_vec_zero(W, Wlen - 1);
    fmpz_one(W + Wlen - 1);
    _fmpz_poly_div_basecase(Qinv, W, W, Wlen, Qrev, Qlen);
    _fmpz_poly_reverse(Qinv, Qinv, n, n);
}

#define MULLOW(z, x, xn, y, yn, nn) \
    if ((xn) >= (yn)) \
        _fmpz_poly_mullow(z, x, xn, y, yn, nn); \
    else \
        _fmpz_poly_mullow(z, y, yn, x, xn, nn); \

void
_fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 1)
    {
        fmpz_set(Qinv, Q);
        _fmpz_vec_zero(Qinv + 1, n - 1);
    }
    else
    {
        slong alloc, Qnlen, Wlen, W2len;
        fmpz * W;

        alloc = FLINT_MAX(n, 3 * FMPZ_POLY_INV_NEWTON_CUTOFF);
        W = _fmpz_vec_init(alloc);

        FLINT_NEWTON_INIT(FMPZ_POLY_INV_NEWTON_CUTOFF, n)

        FLINT_NEWTON_BASECASE(n)
        _fmpz_poly_inv_series_basecase_rev(Qinv, W, Q, Qlen, n);
        FLINT_NEWTON_END_BASECASE

        FLINT_NEWTON_LOOP(m, n)
        Qnlen = FLINT_MIN(Qlen, n);
        Wlen = FLINT_MIN(Qnlen + m - 1, n);
        W2len = Wlen - m;
        MULLOW(W, Q, Qnlen, Qinv, m, Wlen);
        MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m);
        _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);
        FLINT_NEWTON_END_LOOP

        FLINT_NEWTON_END

        _fmpz_vec_clear(W, alloc);
    }
}

void
fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_printf("Exception (fmpz_poly_inv_series_newton). Division by zero.\n");
        flint_abort();
    }

    if (Qinv != Q)
    {
        fmpz_poly_fit_length(Qinv, n);
        _fmpz_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_inv_series_newton(t->coeffs, Q->coeffs, Qlen, n);
        fmpz_poly_swap(Qinv, t);
        fmpz_poly_clear(t);
    }

    _fmpz_poly_set_length(Qinv, n);
    _fmpz_poly_normalise(Qinv);
}

