/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2019 William Hart
    Copyright (C) 2014, 2021, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_div_series_newton(gr_ptr res, gr_srcptr B, slong Blen, gr_srcptr A, slong Alen, slong len, slong cutoff, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, m, n, Anlen, Wlen, W2len, alloc;
    gr_ptr W, T;
    slong a[FLINT_BITS];

    if (len == 0)
        return GR_SUCCESS;

    if (Alen == 0)
        return GR_DOMAIN;

    Alen = FLINT_MIN(Alen, len);
    Blen = FLINT_MIN(Blen, len);

    /* not supported by the following code */
    if (Alen == 1)
        return _gr_poly_div_series_basecase(res, B, Blen, A, Alen, len, ctx);

    cutoff = FLINT_MAX(cutoff, 2);

    a[i = 0] = n = len;
    while (n >= cutoff)
        a[++i] = (n = (n + 1) / 2);

    status |= _gr_poly_inv_series_basecase(res, A, Alen, n, ctx);

    if (status != GR_SUCCESS)
        return status;

    alloc = len + (len + 1) / 2;
    GR_TMP_INIT_VEC(W, alloc, ctx);
    T = GR_ENTRY(W, len, sz);

    for (i--; i >= 1; i--)
    {
        m = n;
        n = a[i];

        Anlen = FLINT_MIN(Alen, n);
        Wlen = FLINT_MIN(Anlen + m - 1, n);
        W2len = Wlen - m;
        status |= _gr_poly_mullow(W, A, Anlen, res, m, Wlen, ctx);
        if (W2len != 0)
            status |= _gr_poly_mullow(GR_ENTRY(res, m, sz), res, m, GR_ENTRY(W, m, sz), W2len, n - m, ctx);
        status |= _gr_vec_neg(GR_ENTRY(res, m, sz), GR_ENTRY(res, m, sz), n - m, ctx);
    }

    m = (len + 1) / 2;
    n = len;

    Anlen = FLINT_MIN(Alen, n);
    Wlen = FLINT_MIN(Anlen + m - 1, n);

    /* Karp-Markstein */
    status |= _gr_poly_mullow(T, res, m, B, Blen, m, ctx);
    status |= _gr_poly_mullow(W, A, Anlen, T, m, Wlen, ctx);
    status |= _gr_poly_sub(GR_ENTRY(W, m, sz), GR_ENTRY(B, m, sz), FLINT_MAX(0, FLINT_MIN(Blen - m, n - m)), GR_ENTRY(W, m, sz), n - m, ctx);
    status |= _gr_poly_mullow(GR_ENTRY(res, m, sz), res, m, GR_ENTRY(W, m, sz), n - m, n - m, ctx);
    _gr_vec_swap(res, T, m, ctx);

    GR_TMP_CLEAR_VEC(W, alloc, ctx);

    return status;
}

int
gr_poly_div_series_newton(gr_poly_t Q, const gr_poly_t A, const gr_poly_t B, slong len, slong cutoff, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (len == 0)
        return gr_poly_zero(Q, ctx);

    if (B->length == 0)
        return GR_DOMAIN;

    if (A->length == 0)
    {
        truth_t is_zero = gr_poly_is_zero(B, ctx);

        if (is_zero == T_FALSE)
            return gr_poly_zero(Q, ctx);

        return GR_UNABLE;
    }

    if (Q == A || Q == B)
    {
        gr_poly_t t;
        gr_poly_init(t, ctx);
        status = gr_poly_div_series_newton(t, A, B, len, cutoff, ctx);
        gr_poly_swap(Q, t, ctx);
        gr_poly_clear(t, ctx);
        return status;
    }

    gr_poly_fit_length(Q, len, ctx);
    status = _gr_poly_div_series_newton(Q->coeffs, A->coeffs, A->length, B->coeffs, B->length, len, cutoff, ctx);
    _gr_poly_set_length(Q, len, ctx);
    _gr_poly_normalise(Q, ctx);
    return status;
}
