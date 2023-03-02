/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#define FQ_POLY_INV_NEWTON_CUTOFF  64

void
_TEMPLATE(T, poly_inv_series_newton) (TEMPLATE(T, struct) * Qinv,
                                      const TEMPLATE(T, struct) * Q, slong n,
                                      const TEMPLATE(T, t) cinv,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    if (n == 1)                 /* {Q,1} x* cinv == 1 mod (x) */
    {
        TEMPLATE(T, set) (Qinv, cinv, ctx);
    }
    else
    {
        const slong alloc = FLINT_MAX(n, 3 * FQ_POLY_INV_NEWTON_CUTOFF);
        slong *a, i, m;
        TEMPLATE(T, struct) * W;

        W = _TEMPLATE(T, vec_init) (alloc, ctx);

        for (i = 1; (WORD(1) << i) < n; i++) ;

        a = (slong *) flint_malloc(i * sizeof(slong));
        a[i = 0] = n;
        while (n >= FQ_POLY_INV_NEWTON_CUTOFF)
            a[++i] = (n = (n + 1) / 2);

        /* Base case */
        {
            TEMPLATE(T, struct) * Qrev = W + 2 * FQ_POLY_INV_NEWTON_CUTOFF;

            _TEMPLATE(T, poly_reverse) (Qrev, Q, n, n, ctx);
            _TEMPLATE(T, vec_zero) (W, 2 * n - 2, ctx);
            TEMPLATE(T, one) (W + (2 * n - 2), ctx);
            _TEMPLATE(T, poly_div_basecase) (Qinv, W, W, 2 * n - 1, Qrev, n,
                                             cinv, ctx);
            _TEMPLATE(T, poly_reverse) (Qinv, Qinv, n, n, ctx);
        }

        for (i--; i >= 0; i--)
        {
            m = n;
            n = a[i];

            _TEMPLATE(T, poly_mullow) (W, Q, n, Qinv, m, n, ctx);
            _TEMPLATE(T, poly_mullow) (Qinv + m, Qinv, m, W + m, n - m, n - m,
                                       ctx);
            _TEMPLATE(T, poly_neg) (Qinv + m, Qinv + m, n - m, ctx);
        }

        _TEMPLATE(T, vec_clear) (W, alloc, ctx);
        flint_free(a);
    }
}

void
TEMPLATE(T, poly_inv_series_newton) (TEMPLATE(T, poly_t) Qinv,
                                     const TEMPLATE(T, poly_t) Q, slong n,
                                     const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, t) cinv;
    TEMPLATE(T, struct) * Qcopy;
    int Qalloc;

    if (Q->length >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        Qcopy = _TEMPLATE(T, vec_init) (n, ctx);
        _TEMPLATE(T, vec_set) (Qcopy, Q->coeffs, Q->length, ctx);
        Qalloc = 1;
    }

    TEMPLATE(T, init) (cinv, ctx);
    TEMPLATE(T, inv) (cinv, Q->coeffs, ctx);

    if (Qinv != Q)
    {
        TEMPLATE(T, poly_fit_length) (Qinv, n, ctx);
        _TEMPLATE(T, poly_inv_series_newton) (Qinv->coeffs, Qcopy, n, cinv,
                                              ctx);
    }
    else
    {
        TEMPLATE(T, struct) * t = _TEMPLATE(T, vec_init) (n, ctx);

        _TEMPLATE(T, poly_inv_series_newton) (t, Qcopy, n, cinv, ctx);

        _TEMPLATE(T, vec_clear) (Qinv->coeffs, Qinv->alloc, ctx);
        Qinv->coeffs = t;
        Qinv->alloc = n;
        Qinv->length = n;
    }
    _TEMPLATE(T, poly_set_length) (Qinv, n, ctx);
    _TEMPLATE(T, poly_normalise) (Qinv, ctx);

    if (Qalloc)
        _TEMPLATE(T, vec_clear) (Qcopy, n, ctx);
    TEMPLATE(T, clear) (cinv, ctx);
}


#endif
