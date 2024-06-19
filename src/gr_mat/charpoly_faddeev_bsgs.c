/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"

static int
gr_mat_trace_prod(gr_ptr res, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong n, i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    n = A->r;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == 0 && j == 0)
                status |= gr_mul(res, GR_MAT_ENTRY(A, 0, 0, sz), GR_MAT_ENTRY(B, 0, 0, sz), ctx);
            else
                status |= gr_addmul(res, GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(B, j, i, sz), ctx);
        }
    }

    return status;
}

static int
gr_mat_trace_prod2(gr_ptr res, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong n, i;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    n = A->r;

    for (i = 0; i < n; i++)
        status |= _gr_vec_dot(res, ((i == 0) ? NULL : res), 0, GR_MAT_ENTRY(A, i, 0, sz), GR_MAT_ENTRY(B, i, 0, sz), n, ctx);

    return status;
}

int
_gr_mat_charpoly_faddeev_bsgs(gr_ptr c, gr_mat_t adj, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t B, C;
    gr_ptr t;
    gr_mat_struct * Apow;
    slong n, m, k, j, i, m_orig;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = A->r;

    if (n == 0)
        return gr_one(c, ctx);

    if (n == 1)
    {
        status |= gr_neg(c, GR_MAT_ENTRY(A, 0, 0, sz), ctx);
        status |= gr_one(GR_ENTRY(c, 1, sz), ctx);

        if (adj != NULL)
            status |= gr_mat_one(adj, ctx);

        return status;
    }

    m_orig = m = n_sqrt(n);

    Apow = flint_malloc(sizeof(gr_mat_struct) * (m + 1));
    t = flint_malloc(sz * (m + 1));
    _gr_vec_init(t, m + 1, ctx);
    for (k = 0; k <= m; k++)
        gr_mat_init(Apow + k, n, n, ctx);

    status |= gr_mat_set(Apow + 1, A, ctx);
    for (k = 2; k <= m; k++)
        status |= gr_mat_mul(Apow + k, Apow + k - 1, Apow + 1, ctx);

    for (k = 1; k <= m; k++)
        status |= gr_mat_trace(GR_ENTRY(t, k, sz), Apow + k, ctx);

    gr_mat_init(B, n, n, ctx);
    gr_mat_init(C, n, n, ctx);

    status |= gr_one(GR_ENTRY(c, n, sz), ctx);
    status |= gr_mat_trace(GR_ENTRY(c, n - 1, sz), A, ctx);
    status |= gr_neg(GR_ENTRY(c, n - 1, sz), GR_ENTRY(c, n - 1, sz), ctx);
    status |= gr_mat_add_scalar(B, A, GR_ENTRY(c, n - 1, sz), ctx);

    k = 2;
    while (k <= n - 1)
    {
        m = FLINT_MIN(m, n - k);

        status |= gr_mat_transpose(B, B, ctx);

        status |= gr_mat_trace_prod2(GR_ENTRY(c, n - k, sz), A, B, ctx);
        status |= gr_div_si(GR_ENTRY(c, n - k, sz), GR_ENTRY(c, n - k, sz), -k, ctx);
        if (status != GR_SUCCESS)
            goto cleanup;

        for (j = 1; j < m; j++)
        {
            status |= gr_mat_trace_prod2(GR_ENTRY(c, n - k - j, sz), Apow + j + 1, B, ctx);

            for (i = 0; i < j; i++)
                status |= gr_addmul(GR_ENTRY(c, n - k - j, sz), GR_ENTRY(t, j - i, sz), GR_ENTRY(c, n - k - i, sz), ctx);

            status |= gr_div_si(GR_ENTRY(c, n - k - j, sz), GR_ENTRY(c, n - k - j, sz), -k - j, ctx);
            if (status != GR_SUCCESS)
                goto cleanup;
        }

        status |= gr_mat_transpose(B, B, ctx);

        status |= gr_mat_mul(C, Apow + m, B, ctx);
        gr_mat_swap(B, C, ctx);

        for (j = 0; j < m; j++)
        {
            if (m - j - 1 == 0)
                status |= gr_mat_add_scalar(B, B, GR_ENTRY(c, n - k - j, sz), ctx);
            else
                status |= gr_mat_addmul_scalar(B, Apow + m - j - 1, GR_ENTRY(c, n - k - j, sz), ctx);
        }

        k += m;
    }

    status |= gr_mat_trace_prod(c, A, B, ctx);
    status |= gr_div_si(c, c, -n, ctx);

    if (adj != NULL)
    {
        if (n % 2)
            status |= gr_mat_set(adj, B, ctx);
        else
            status |= gr_mat_neg(adj, B, ctx);
    }

cleanup:
    for (k = 0; k <= m_orig; k++)
        gr_mat_clear(Apow + k, ctx);
    flint_free(Apow);
    _gr_vec_clear(t, m_orig + 1, ctx);
    flint_free(t);

    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
gr_mat_charpoly_faddeev_bsgs(gr_poly_t cp, gr_mat_t adjugate, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_faddeev_bsgs(cp->coeffs, adjugate, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
