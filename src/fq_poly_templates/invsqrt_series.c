/*
    Copyright (C) 2011, 2021, 2022 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

static void
_TEMPLATE(T, poly_invsqrt_series_prealloc)(TEMPLATE(T, struct) * g,
                        const TEMPLATE(T, struct) * h, TEMPLATE(T, struct) * t,
                        TEMPLATE(T, struct) * u, slong n, const TEMPLATE(T, ctx_t) ctx)
{
    const int alloc = (t == NULL);
    const slong m   = (n + 1) / 2;
    TEMPLATE(T, t) c, inv2, one;

    if (n == 1)
    {
        TEMPLATE(T, set_ui)(g + 0, 1, ctx);
        return;
    }

    if (alloc)
    {
        t = _TEMPLATE(T, vec_init)(n, ctx);
        u = _TEMPLATE(T, vec_init)(n, ctx);
    }

    TEMPLATE(T, init)(c, ctx);
    TEMPLATE(T, init)(inv2, ctx);
    TEMPLATE(T, init)(one, ctx);

    TEMPLATE(T, set_ui)(one, 1, ctx);

    TEMPLATE(T, set_ui)(inv2, 2, ctx);
    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
       TEMPLATE(T, inv)(inv2, inv2, ctx);

    _TEMPLATE(T, poly_invsqrt_series_prealloc)(g, h, t, u, m, ctx);

    _TEMPLATE(T, vec_zero)(g + m, n - m, ctx);

    _TEMPLATE(T, poly_mul)(t, g, m, g, m, ctx);
    if (2*m - 1 < n)
        TEMPLATE(T, zero)(t + n - 1, ctx);

    _TEMPLATE(T, poly_mullow)(u, t, n, g, n, n, ctx);
    _TEMPLATE(T, poly_mullow)(t, u, n, h, n, n, ctx);

    TEMPLATE(T, sub)(c, c, one, ctx);
    TEMPLATE(T, mul)(c, c, inv2, ctx);
    _TEMPLATE3(T, vec_scalar_mul, T)(g + m, t + m, n - m, c, ctx);

    if (alloc)
    {
        _TEMPLATE(T, vec_clear)(t, n, ctx);
        _TEMPLATE(T, vec_clear)(u, n, ctx);
    }

    TEMPLATE(T, clear)(one, ctx);
    TEMPLATE(T, clear)(inv2, ctx);
    TEMPLATE(T, clear)(c, ctx);
}

void _TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, struct) * g,
           const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    _TEMPLATE(T, poly_invsqrt_series_prealloc)(g, h, NULL, NULL, n, ctx);
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
        flint_throw(FLINT_ERROR, "Exception (fq_poly_invsqrt). Division by zero.\n");
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

