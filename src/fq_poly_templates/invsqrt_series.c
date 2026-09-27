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

/*
    Third-order Newton iteration: given g = h^(-1/2) + O(x^m), let
    e = 1 - h g^2 = O(x^m). Then

        h^(-1/2) = g (1 + e/2 + 3e^2/8) + O(x^(3m)).

    Assumes that h has length n.
*/
void _TEMPLATE(T, poly_invsqrt_series)(TEMPLATE(T, struct) * g,
           const TEMPLATE(T, struct) * h, slong n, TEMPLATE(T, ctx_t) ctx)
{
    slong a[FLINT_BITS];
    slong i, m, k, L, len, tlen;
    TEMPLATE(T, struct) * t, * u, * v, * w;
    TEMPLATE(T, t) inv2, inv8;

    TEMPLATE(T, one)(g, ctx);

    if (n == 1)
        return;

    len = n;
    a[i = 0] = n;
    while (n > 1)
        a[++i] = (n = (n + 2) / 3);

    t = _TEMPLATE(T, vec_init)(3 * len, ctx);
    u = t + len;
    v = u + len;

    TEMPLATE(T, init)(inv2, ctx);
    TEMPLATE(T, init)(inv8, ctx);

    /* -1/2 and -1/8; these are zero (the result is undefined) in
       characteristic 2 */
#if defined(FQ_NMOD_POLY_H) || defined(FQ_ZECH_POLY_H)
    if (TEMPLATE(T, ctx_prime)(ctx) != 2)
#else
    if (fmpz_cmp_ui(TEMPLATE(T, ctx_prime)(ctx), 2) != 0)
#endif
    {
        TEMPLATE(T, set_ui)(inv2, 2, ctx);
        TEMPLATE(T, inv)(inv2, inv2, ctx);
        TEMPLATE(T, neg)(inv2, inv2, ctx);
        TEMPLATE(T, set_ui)(inv8, 8, ctx);
        TEMPLATE(T, inv)(inv8, inv8, ctx);
        TEMPLATE(T, neg)(inv8, inv8, ctx);
    }

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];
        k = n - 2 * m;
        L = n - m;

        /* w = (h g^2)[m, n) = -e / x^m */
        tlen = FLINT_MIN(2 * m - 1, n);
        _TEMPLATE(T, poly_mullow)(t, g, m, g, m, tlen, ctx);
        _TEMPLATE(T, poly_mullow)(u, t, tlen, h, n, n, ctx);
        w = u + m;

        if (k > 0)
        {
            /* e/2 + 3e^2/8 = -(4w - 3 x^m w^2) / 8 */
            _TEMPLATE(T, poly_mullow)(v, w, k, w, k, k, ctx);
            _TEMPLATE(T, vec_scalar_mul_ui)(v, v, k, 3, ctx);
            _TEMPLATE(T, vec_scalar_mul_ui)(w, w, L, 4, ctx);
            _TEMPLATE(T, vec_sub)(w + m, w + m, v, k, ctx);
        }

        _TEMPLATE(T, poly_mullow)(g + m, g, FLINT_MIN(m, L), w, L, L, ctx);
        _TEMPLATE3(T, vec_scalar_mul, T)(g + m, g + m, L, (k > 0) ? inv8 : inv2, ctx);
    }

    TEMPLATE(T, clear)(inv2, ctx);
    TEMPLATE(T, clear)(inv8, ctx);
    _TEMPLATE(T, vec_clear)(t, 3 * len, ctx);
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
