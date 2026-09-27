/*
    Copyright (C) 2011, 2021, 2022 William Hart
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

/*
    Karp-Markstein: with g = h^(-1/2) + O(x^m) and v = h g mod x^m,
    sqrt(h) = v + g (h - v^2) / 2 + O(x^(2m)).

    Assumes that h has length n.
*/
void
_TEMPLATE(T, poly_sqrt_series)(TEMPLATE(T, struct) * g, const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * r, * t;
    TEMPLATE(T, t) inv2;
    slong m, L, tlen;

    if (n <= 2)
    {
        t = _TEMPLATE(T, vec_init)(n, ctx);
        _TEMPLATE(T, poly_invsqrt_series)(t, h, n, ctx);
        _TEMPLATE(T, poly_mullow)(g, t, n, h, n, n, ctx);
        _TEMPLATE(T, vec_clear)(t, n, ctx);
        return;
    }

    m = (n + 1) / 2;
    L = n - m;
    tlen = FLINT_MIN(2 * m - 1, n);

    r = _TEMPLATE(T, vec_init)(m + n, ctx);
    t = r + m;

    TEMPLATE(T, init)(inv2, ctx);
#if defined(FQ_NMOD_POLY_H) || defined(FQ_ZECH_POLY_H)
    if (TEMPLATE(T, ctx_prime)(ctx) != 2)
#else
    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
#endif
    {
        TEMPLATE(T, set_ui)(inv2, 2, ctx);
        TEMPLATE(T, inv)(inv2, inv2, ctx);
    }

    _TEMPLATE(T, poly_invsqrt_series)(r, h, m, ctx);
    _TEMPLATE(T, poly_mullow)(g, r, m, h, m, m, ctx);

    /* t[m, n) = h[m, n) - (v^2)[m, n) */
    _TEMPLATE(T, poly_mullow)(t, g, m, g, m, tlen, ctx);
    _TEMPLATE(T, vec_sub)(t + m, h + m, t + m, tlen - m, ctx);
    _TEMPLATE(T, vec_set)(t + tlen, h + tlen, n - tlen, ctx);

    _TEMPLATE(T, poly_mullow)(g + m, r, L, t + m, L, L, ctx);
    _TEMPLATE3(T, vec_scalar_mul, T)(g + m, g + m, L, inv2, ctx);

    TEMPLATE(T, clear)(inv2, ctx);
    _TEMPLATE(T, vec_clear)(r, m + n, ctx);
}

void
TEMPLATE(T, poly_sqrt_series)(TEMPLATE(T, poly_t) g, const TEMPLATE(T, poly_t) h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, struct) * g_coeffs;
    TEMPLATE(T, struct) * h_coeffs;
    TEMPLATE(T, poly_t) t1;
    slong hlen;

    hlen = h->length;

    if (n == 0)
    {
        flint_throw(FLINT_DIVZERO, "Exception (fq_poly_sqrt_series). Division by zero.\n");
    }

    if (h->length == 0 || !TEMPLATE(T, is_one)(h->coeffs + 0, ctx))
    {
        flint_throw(FLINT_ERROR, "Exception (fq_poly_sqrt_series). Requires constant term 1.\n");
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

    _TEMPLATE(T, poly_sqrt_series)(g_coeffs, h_coeffs, n, ctx);

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
