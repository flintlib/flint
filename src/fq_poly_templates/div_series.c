/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2014 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"
#include "gr_poly.h"

/* Below this length, the basecase algorithm is faster than Newton
   iteration. */
#if defined(FQ_ZECH_POLY_H)
#define FQ_POLY_DIV_SERIES_NEWTON_CUTOFF 64
#elif defined(FQ_NMOD_POLY_H)
#define FQ_POLY_DIV_SERIES_NEWTON_CUTOFF 16
#else
#define FQ_POLY_DIV_SERIES_NEWTON_CUTOFF 32
#endif

void
_TEMPLATE(T, poly_div_series) (TEMPLATE(T, struct) * Q, const TEMPLATE(T, struct) * A, slong Alen,
    const TEMPLATE(T, struct) * B, slong Blen, slong n, const TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;
    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);

    if (n < FQ_POLY_DIV_SERIES_NEWTON_CUTOFF || FLINT_MIN(Blen, n) < 10)
    {
        TEMPLATE(T, t) u;
        TEMPLATE(T, init)(u, ctx);
        TEMPLATE(T, inv)(u, B + 0, ctx);
        GR_MUST_SUCCEED(_gr_poly_div_series_basecase_preinv1(Q, A, Alen, B, Blen, u, n, gr_ctx));
        TEMPLATE(T, clear)(u, ctx);
    }
    else
    {
        if (TEMPLATE(T, is_zero)(B + 0, ctx))
            TEMPLATE(T, inv)(Q, B + 0, ctx);  /* throws */

        GR_MUST_SUCCEED(_gr_poly_div_series_newton(Q, A, Alen, B, Blen, n, FQ_POLY_DIV_SERIES_NEWTON_CUTOFF, gr_ctx));
    }
}

void TEMPLATE(T, poly_div_series)(TEMPLATE(T, poly_t) Q, const TEMPLATE(T, poly_t) A,
                                         const TEMPLATE(T, poly_t) B, slong n, const TEMPLATE(T, ctx_t) ctx)
{
    slong Alen = FLINT_MIN(A->length, n);
    slong Blen = FLINT_MIN(B->length, n);

    if (Blen == 0)
    {
        flint_throw(FLINT_DIVZERO, "Exception (fq_poly_div_series). Division by zero.\n");
    }

    if (Alen == 0)
    {
        TEMPLATE(T, poly_zero)(Q, ctx);
        return;
    }

    if (Q == A || Q == B)
    {
        TEMPLATE(T, poly_t) t;
        TEMPLATE(T, poly_init2)(t, n, ctx);
        _TEMPLATE(T, poly_div_series)(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, ctx);
        TEMPLATE(T, poly_swap)(Q, t, ctx);
        TEMPLATE(T, poly_clear)(t, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length)(Q, n, ctx);
        _TEMPLATE(T, poly_div_series)(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, n, ctx);
    }

    _TEMPLATE(T, poly_set_length)(Q, n, ctx);
    _TEMPLATE(T, poly_normalise)(Q, ctx);
}

#endif
