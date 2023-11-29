/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define MULLOW(z, x, xn, y, yn, nn) \
    if ((xn) >= (yn)) \
        _fmpz_poly_mullow(z, x, xn, y, yn, nn); \
    else \
        _fmpz_poly_mullow(z, y, yn, x, xn, nn); \

void
_fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
{
    slong cutoff = 64;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen < cutoff)
    {
        _fmpz_poly_inv_series_basecase(Qinv, Q, Qlen, n);
    }
    else
    {
        slong *a, i, m, Qnlen, Wlen, W2len;
        fmpz * W;

        W = _fmpz_vec_init(n);
        a = flint_malloc(sizeof(slong) * FLINT_BITS);

        a[i = 0] = n;
        while (n >= cutoff)
            a[++i] = (n = (n + 1) / 2);

        _fmpz_poly_inv_series_basecase(Qinv, Q, Qlen, n);

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            Qnlen = FLINT_MIN(Qlen, n);
            Wlen = FLINT_MIN(Qnlen + m - 1, n);
            W2len = Wlen - m;
            MULLOW(W, Q, Qnlen, Qinv, m, Wlen);
            MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m);
            _fmpz_vec_neg(Qinv + m, Qinv + m, n - m);
        }

        _fmpz_vec_clear(W, n);
        flint_free(a);
    }
}

void
fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_inv_series_newton). Division by zero.\n");
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

