/*
    Copyright (C) 2011, 2021, 2022 William Hart
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2024 Albin Ahlbäck

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
#define FQ_POLY_INVSQRT_SERIES_NEWTON_CUTOFF 64
#elif defined(FQ_NMOD_POLY_H)
#define FQ_POLY_INVSQRT_SERIES_NEWTON_CUTOFF 12
#else
#define FQ_POLY_INVSQRT_SERIES_NEWTON_CUTOFF 16
#endif

void _TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, struct) * g,
           const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    gr_ctx_t gr_ctx;

#if defined(FQ_NMOD_POLY_H) || defined(FQ_ZECH_POLY_H)
    if (n > 1 && TEMPLATE(T, ctx_prime)(ctx) == 2)
#else
    if (n > 1 && fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) == 0)
#endif
        flint_throw(FLINT_ERROR, "(%s): Not supported in characteristic 2.\n", __func__);

    TEMPLATE3(_gr_ctx_init, T, from_ref)(gr_ctx, ctx);
    GR_MUST_SUCCEED(_gr_poly_rsqrt_series_newton(g, h, n, n, FQ_POLY_INVSQRT_SERIES_NEWTON_CUTOFF, gr_ctx));
}

void TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, poly_t) g,
               const TEMPLATE(T, poly_t) h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    const slong hlen = h->length;
    TEMPLATE(T, struct) * g_coeffs;
    TEMPLATE(T, struct) * h_coeffs;
    TEMPLATE(T, poly_t) t1;

    if (n == 0 || h->length == 0 || TEMPLATE(T, is_zero)(h->coeffs + 0, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Exception (fq_poly_invsqrt). Division by zero.\n");
    }

    if (!TEMPLATE(T, is_one)(h->coeffs + 0, ctx))
    {
        flint_throw(FLINT_ERROR, "Exception (fq_poly_invsqrt_series). Constant term != 1.\n");
    }

    if (hlen < n)
    {
        h_coeffs = _TEMPLATE(T, vec_init)(n, ctx);
        _TEMPLATE(T, vec_set)(h_coeffs, h->coeffs, hlen, ctx);
    }
    else
        h_coeffs = h->coeffs;

    if (h == g && hlen >= n)
    {
        TEMPLATE(T, poly_init2)(t1, n, ctx);
        g_coeffs = t1->coeffs;
    }
    else
    {
        TEMPLATE(T, poly_fit_length)(g, n, ctx);
        g_coeffs = g->coeffs;
    }

    _TEMPLATE(T, poly_invsqrt_series)(g_coeffs, h_coeffs, n, ctx);

    if (h == g && hlen >= n)
    {
        TEMPLATE(T, poly_swap)(g, t1, ctx);
        TEMPLATE(T, poly_clear)(t1, ctx);
    }

    g->length = n;

    if (hlen < n)
        _TEMPLATE(T, vec_clear)(h_coeffs, n, ctx);

    _TEMPLATE(T, poly_normalise)(g, ctx);
}

#endif
