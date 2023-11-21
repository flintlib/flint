/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

#define MULLOW(z, x, xn, y, yn, nn, prec) \
    if ((xn) >= (yn)) \
        _arb_poly_mullow(z, x, xn, y, yn, nn, prec); \
    else \
        _arb_poly_mullow(z, y, yn, x, xn, nn, prec); \

void
_arb_poly_inv_series(arb_ptr Qinv,
    arb_srcptr Q, slong Qlen, slong len, slong prec)
{
    Qlen = FLINT_MIN(Qlen, len);

    arb_inv(Qinv, Q, prec);

    if (Qlen == 1)
    {
        _arb_vec_zero(Qinv + 1, len - 1);
    }
    else if (len == 2)
    {
        arb_mul(Qinv + 1, Qinv, Qinv, prec);
        arb_mul(Qinv + 1, Qinv + 1, Q + 1, prec);
        arb_neg(Qinv + 1, Qinv + 1);
    }
    else
    {
        slong i, blen;

        /* The basecase algorithm is faster for much larger Qlen or len than
           this, but unfortunately also much less numerically stable. */
        if (Qlen == 2 || len <= 8)
            blen = len;
        else
            blen = FLINT_MIN(len, 4);

        for (i = 1; i < blen; i++)
        {
            arb_dot(Qinv + i, NULL, 1,
                Q + 1, 1, Qinv + i - 1, -1, FLINT_MIN(i, Qlen - 1), prec);
            if (!arb_is_one(Qinv))
                arb_mul(Qinv + i, Qinv + i, Qinv, prec);
        }

        if (len > blen)
        {
            slong Qnlen, Wlen, W2len;
            arb_ptr W;

            W = _arb_vec_init(len);

            NEWTON_INIT(blen, len)
            NEWTON_LOOP(m, n)

            Qnlen = FLINT_MIN(Qlen, n);
            Wlen = FLINT_MIN(Qnlen + m - 1, n);
            W2len = Wlen - m;
            MULLOW(W, Q, Qnlen, Qinv, m, Wlen, prec);
            MULLOW(Qinv + m, Qinv, m, W + m, W2len, n - m, prec);
            _arb_vec_neg(Qinv + m, Qinv + m, n - m);

            NEWTON_END_LOOP
            NEWTON_END

            _arb_vec_clear(W, len);
        }
    }
}

void
arb_poly_inv_series(arb_poly_t Qinv, const arb_poly_t Q, slong n, slong prec)
{
    if (n == 0)
    {
        arb_poly_zero(Qinv);
        return;
    }

    if (Q->length == 0)
    {
        arb_poly_fit_length(Qinv, n);
        _arb_vec_indeterminate(Qinv->coeffs, n);
        _arb_poly_set_length(Qinv, n);
        return;
    }

    if (Qinv == Q)
    {
        arb_poly_t t;
        arb_poly_init(t);
        arb_poly_inv_series(t, Q, n, prec);
        arb_poly_swap(Qinv, t);
        arb_poly_clear(t);
        return;
    }

    arb_poly_fit_length(Qinv, n);
    _arb_poly_inv_series(Qinv->coeffs, Q->coeffs, Q->length, n, prec);
    _arb_poly_set_length(Qinv, n);
    _arb_poly_normalise(Qinv);
}

