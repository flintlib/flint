/*
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

/* Below this length, the basecase algorithm (using the precomputed
   inverse of the constant term) is faster than Newton iteration. */
#if defined(FQ_ZECH_POLY_H)
#define FQ_POLY_INV_NEWTON_CUTOFF 64
#elif defined(FQ_NMOD_POLY_H)
#define FQ_POLY_INV_NEWTON_CUTOFF 16
#else
#define FQ_POLY_INV_NEWTON_CUTOFF 32
#endif

void
_TEMPLATE(T, poly_inv_series_newton) (TEMPLATE(T, struct) * Qinv,
                                      const TEMPLATE(T, struct) * Q, slong n,
                                      const TEMPLATE(T, t) cinv,
                                      const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    if (n < FQ_POLY_INV_NEWTON_CUTOFF)
        GR_MUST_SUCCEED(_gr_poly_inv_series_basecase_preinv1(Qinv, Q, n, cinv, n, gr_ctx));
    else
        GR_MUST_SUCCEED(_gr_poly_inv_series_newton(Qinv, Q, n, n, FQ_POLY_INV_NEWTON_CUTOFF, gr_ctx));
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
