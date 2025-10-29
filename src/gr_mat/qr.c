/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_lq_gso(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong m, n, i, k;
    gr_ptr Qi, Qk, s, t;

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);

    if (m > n || m != gr_mat_nrows(L, ctx) || m != gr_mat_ncols(L, ctx)
              || m != gr_mat_nrows(Q, ctx) || n != gr_mat_ncols(Q, ctx))
    {
        return GR_DOMAIN;
    }

    if (A == L)
    {
        gr_mat_t T;
        gr_mat_init(T, m, m, ctx);
        status = gr_mat_lq_gso(T, Q, A, ctx);
        status |= gr_mat_swap_entrywise(L, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    GR_TMP_INIT(t, ctx);

    for (k = 0; k < m && status == GR_SUCCESS; k++)
    {
        Qk = GR_MAT_ENTRY(Q, k, 0, sz);

        if (A != Q)
            status |= _gr_vec_set(Qk, GR_MAT_ENTRY(A, k, 0, sz), n, ctx);

        for (i = 0; i < k; i++)
        {
            Qi =  GR_MAT_ENTRY(Q, i, 0, sz);
            s = GR_MAT_ENTRY(L, k, i, sz);
            status |= _gr_vec_dot(s, NULL, 0, Qi, Qk, n, ctx);
            status |= _gr_vec_submul_scalar(Qk, Qi, n, s, ctx);
        }

        s = GR_MAT_ENTRY(L, k, k, sz);
        /* todo: optimized vec_norm2 (or optimize self-dot for nfloat, arf etc) */
        status |= _gr_vec_dot(s, NULL, 0, Qk, Qk, n, ctx);
        status |= gr_sqrt(s, s, ctx);
        status |= gr_inv(t, s, ctx);
        status |= _gr_vec_mul_scalar(Qk, Qk, n, t, ctx);
        status |= _gr_vec_zero(GR_MAT_ENTRY(L, k, k + 1, sz), m - k - 1, ctx);
    }

    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_mat_lq_recursive(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong m, n, d;
    gr_mat_t A1, A2, Q1, Q2, L11, L12, L21, L22, QT, T;

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);

    if (m <= 1)
        return gr_mat_lq_gso(L, Q, A, ctx);

    if (m > n || m != gr_mat_nrows(L, ctx) || m != gr_mat_ncols(L, ctx)
              || m != gr_mat_nrows(Q, ctx) || n != gr_mat_ncols(Q, ctx))
    {
        return GR_DOMAIN;
    }

    if (A == L)
    {
        gr_mat_t T;
        gr_mat_init(T, m, m, ctx);
        status = gr_mat_lq_recursive(T, Q, A, ctx);
        status |= gr_mat_swap_entrywise(L, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    d = m / 2;

    gr_mat_window_init(A1, A, 0, 0, d, n, ctx);
    gr_mat_window_init(A2, A, d, 0, m, n, ctx);
    gr_mat_window_init(Q1, Q, 0, 0, d, n, ctx);
    gr_mat_window_init(Q2, Q, d, 0, m, n, ctx);

    gr_mat_window_init(L11, L, 0, 0, d, d, ctx);
    gr_mat_window_init(L12, L, 0, d, d, m, ctx);
    gr_mat_window_init(L21, L, d, 0, m, d, ctx);
    gr_mat_window_init(L22, L, d, d, m, m, ctx);

    /* A1 = L11 Q1 */
    if (A == Q)
        status |= gr_mat_lq(L11, Q1, Q1, ctx);
    else
        status |= gr_mat_lq(L11, Q1, A1, ctx);

    if (status != GR_SUCCESS)
        return status;

    status |= gr_mat_zero(L12, ctx);

    /* L21 = A2 Q1^T -- todo: mul_transpose */
    GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(QT, Q1, ctx);
    status |= gr_mat_mul(L21, A2, QT, ctx);
    GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(QT, ctx);

    /* Q2 = A2 - L21 Q1 -- todo: in-place submul */
    gr_mat_init(T, Q2->r, Q2->c, ctx);
    status |= gr_mat_mul(T, L21, Q1, ctx);
    status |= gr_mat_sub(Q2, A2, T, ctx);
    gr_mat_clear(T, ctx);

    /* Q2 = L22 Q2 */
    status |= gr_mat_lq(L22, Q2, Q2, ctx);

    return status;
}

int
gr_mat_lq_generic(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
{
    if (A->r <= 4)
        return gr_mat_lq_gso(L, Q, A, ctx);
    else
        return gr_mat_lq_recursive(L, Q, A, ctx);
}

int
gr_mat_lq(gr_mat_t L, gr_mat_t Q, const gr_mat_t A, gr_ctx_t ctx)
{
    return GR_MAT_BINARY_UNARY_OP(ctx, MAT_LQ)(L, Q, A, ctx);
}

int
gr_mat_qr(gr_mat_t Q, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t QT, AT;
    slong m, n;
    int status;

    n = gr_mat_nrows(A, ctx);
    m = gr_mat_ncols(A, ctx);

    if (m > n || m != gr_mat_nrows(R, ctx) || m != gr_mat_ncols(R, ctx)
              || n != gr_mat_nrows(Q, ctx) || m != gr_mat_ncols(Q, ctx))
    {
        return GR_DOMAIN;
    }

    if (A == R)
    {
        gr_mat_t T;
        gr_mat_init(T, m, m, ctx);
        status = gr_mat_qr(Q, T, A, ctx);
        status |= gr_mat_swap_entrywise(R, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    if (A == Q)
    {
        GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(QT, Q, ctx);

        status = gr_mat_lq_gso(R, QT, QT, ctx);

        status |= gr_mat_transpose(R, R, ctx);
        GR_MAT_SHALLOW_TRANSPOSE(Q, QT, ctx);
        GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(QT, ctx);
    }
    else
    {
        GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(AT, A, ctx);
        GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(QT, Q, ctx);

        status = gr_mat_lq_gso(R, QT, AT, ctx);

        status |= gr_mat_transpose(R, R, ctx);
        GR_MAT_SHALLOW_TRANSPOSE(Q, QT, ctx);
        GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(QT, ctx);
        GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(AT, ctx);
    }

    return status;
}

