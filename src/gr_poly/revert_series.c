/*
    Copyright (C) 2011, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"
#include "ulong_extras.h"

#define SPECIAL_CASES \
    if (flen <= 1) \
        return GR_DOMAIN; \
    is_zero = gr_is_zero(f, ctx); \
    if (is_zero == T_UNKNOWN) \
        return GR_UNABLE; \
    if (is_zero == T_FALSE) \
        return GR_DOMAIN; \
    if (n <= 2) \
    { \
        if (n >= 1) \
            status |= gr_zero(res, ctx); \
        if (n == 2) \
            status |= gr_inv(GR_ENTRY(res, 1, sz), GR_ENTRY(f, 1, sz), ctx); \
        return status; \
    }

int
_gr_poly_revert_series_lagrange(gr_ptr res,
    gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
{
    slong i;
    gr_ptr R, S, T;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    truth_t is_zero;

    SPECIAL_CASES

    GR_TMP_INIT_VEC(R, 3 * (n - 1), ctx);
    S = GR_ENTRY(R, n - 1, sz);
    T = GR_ENTRY(S, n - 1, sz);

    status |= gr_zero(res, ctx);
    status |= _gr_poly_inv_series(R, GR_ENTRY(f, 1, sz), FLINT_MIN(flen, n) - 1, n - 1, ctx);
    /* Inverting leading coefficient failed => domain error */
    if (status != GR_SUCCESS)
        goto cleanup;

    status |= gr_set(GR_ENTRY(res, 1, sz), GR_ENTRY(R, 0, sz), ctx);
    status |= _gr_vec_set(S, R, n - 1, ctx);

    for (i = 2; i < n; i++)
    {
        status |= _gr_poly_mullow(T, S, n - 1, R, n - 1, n - 1, ctx);
        status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(T, i - 1, sz), i, ctx);
        FLINT_SWAP(gr_ptr, S, T);
    }

    /* Dividing by i failed in finite characteristic */
    if (status != GR_SUCCESS)
        status = GR_UNABLE;

cleanup:
    GR_TMP_CLEAR_VEC(R, 3 * (n - 1), ctx);
    return status;
}

/* pointer to (x/f)^i */
#define Ri(ii) GR_ENTRY(R, (n-1)*((ii)-1), sz)

int
_gr_poly_revert_series_lagrange_fast(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
{
    slong i, j, m;
    gr_ptr R, S, T, t;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    truth_t is_zero;

    SPECIAL_CASES

    m = n_sqrt(n);

    GR_TMP_INIT_VEC(R, (n - 1) * (m + 2) + 1, ctx);
    S = GR_ENTRY(R, (n - 1) * m, sz);
    T = GR_ENTRY(S, n - 1, sz);
    t = GR_ENTRY(T, n - 1, sz);

    status |= gr_zero(res, ctx);
    status |= _gr_poly_inv_series(Ri(1), GR_ENTRY(f, 1, sz), FLINT_MIN(flen, n) - 1, n - 1, ctx);
    /* Inverting leading coefficient failed => domain error */
    if (status != GR_SUCCESS)
        goto cleanup;

    /* TODO: when do we prefer squaring for powers? (see also compose_series_brent_kung) */
    for (i = 2; i <= m; i++)
        status |= _gr_poly_mullow(Ri(i), Ri((i + 1) / 2), n - 1, Ri(i / 2), n - 1, n - 1, ctx);

    status |= gr_set(GR_ENTRY(res, 1, sz), Ri(1), ctx);
    for (i = 2; i < m; i++)
        status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(Ri(i), i - 1, sz), i, ctx);

    status |= _gr_vec_set(S, Ri(m), n - 1, ctx);

    for (i = m; i < n; i += m)
    {
        status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(S, i - 1, sz), i, ctx);

        for (j = 1; j < m && i + j < n; j++)
        {
            status |= _gr_vec_dot_rev(t, NULL, 0, S, Ri(j), i + j, ctx);
            status |= gr_div_ui(GR_ENTRY(res, i + j, sz), t, i + j, ctx);
        }

        if (i + 1 < n)
        {
            status |= _gr_poly_mullow(T, S, n - 1, Ri(m), n - 1, n - 1, ctx);
            FLINT_SWAP(gr_ptr, S, T);
        }
    }

    /* Dividing by i failed in finite characteristic */
    if (status != GR_SUCCESS)
        status = GR_UNABLE;

cleanup:
    GR_TMP_CLEAR_VEC(R, (n - 1) * (m + 2) + 1, ctx);
    return status;
}

/* Currently, we choose this so that the basecase succeeds in
   small characteristic. */
#define NEWTON_CUTOFF 2

int
_gr_poly_revert_series_newton(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
{
    slong i, k, a[FLINT_BITS];
    gr_ptr T, U, V;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    truth_t is_zero;

    SPECIAL_CASES

    GR_TMP_INIT_VEC(T, 3 * n, ctx);
    U = GR_ENTRY(T, n, sz);
    V = GR_ENTRY(U, n, sz);

    k = n;
    for (i = 1; (WORD(1) << i) < k; i++);
    a[i = 0] = k;
    while (k > NEWTON_CUTOFF)
        a[++i] = (k = (k + 1) / 2);

    status |= _gr_poly_revert_series_lagrange(res, f, flen, k, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    status |= _gr_vec_zero(GR_ENTRY(res, k, sz), n - k, ctx);

    for (i--; i >= 0; i--)
    {
        k = a[i];
        status |= _gr_poly_compose_series(T, f, FLINT_MIN(flen, k), res, k, k, ctx);
        status |= _gr_poly_derivative(U, T, k, ctx);
        status |= gr_zero(GR_ENTRY(U, k - 1, sz), ctx);
        status |= gr_zero(GR_ENTRY(T, 1, sz), ctx);
        status |= _gr_poly_div_series(V, T, k, U, k, k, ctx);
        status |= _gr_poly_derivative(T, res, k, ctx);
        status |= _gr_poly_mullow(U, V, k, T, k, k, ctx);
        status |= _gr_vec_sub(res, res, U, k, ctx);
    }

cleanup:
    GR_TMP_CLEAR_VEC(T, 3 * n, ctx);

    return status;
}

int
_gr_poly_revert_series(gr_ptr res, gr_srcptr f, slong flen, slong n, gr_ctx_t ctx)
{
    int status;

    status = _gr_poly_revert_series_lagrange_fast(res, f, flen, n, ctx);

    /* Newton may succeed where Lagrange fails in finite characteristic */
    if (status == GR_UNABLE)
        status = _gr_poly_revert_series_newton(res, f, flen, n, ctx);

    return status;
}

int
gr_poly_revert_series_lagrange(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (res != f)
    {
        gr_poly_fit_length(res, n, ctx);
        status |= _gr_poly_revert_series_lagrange(res->coeffs, f->coeffs, flen, n, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status |= _gr_poly_revert_series_lagrange(t->coeffs, f->coeffs, flen, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_revert_series_lagrange_fast(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (res != f)
    {
        gr_poly_fit_length(res, n, ctx);
        status |= _gr_poly_revert_series_lagrange_fast(res->coeffs, f->coeffs, flen, n, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status |= _gr_poly_revert_series_lagrange_fast(t->coeffs, f->coeffs, flen, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_revert_series_newton(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (res != f)
    {
        gr_poly_fit_length(res, n, ctx);
        status |= _gr_poly_revert_series_newton(res->coeffs, f->coeffs, flen, n, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status |= _gr_poly_revert_series_newton(t->coeffs, f->coeffs, flen, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}

int
gr_poly_revert_series(gr_poly_t res, const gr_poly_t f, slong n, gr_ctx_t ctx)
{
    slong flen = f->length;
    int status = GR_SUCCESS;

    if (res != f)
    {
        gr_poly_fit_length(res, n, ctx);
        status |= _gr_poly_revert_series(res->coeffs, f->coeffs, flen, n, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status |= _gr_poly_revert_series(t->coeffs, f->coeffs, flen, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
