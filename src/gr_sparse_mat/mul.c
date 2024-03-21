/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2020 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "gr_sparse_mat.h"

int gr_csr_mat_mul_vec(gr_ptr v, const gr_csr_mat_t A, gr_srcptr u, gr_ctx_t ctx)
{
    slong ar, i, sz;
    int status;
    gr_sparse_vec_t row;

    sz = ctx->sizeof_elem;
    ar = gr_sparse_mat_nrows(A, ctx);

    if (gr_csr_mat_is_zero(A, ctx) == T_TRUE)
    {
        return _gr_vec_zero(v, ar, ctx);
    }

    if (u == v)
    {
        gr_ptr w;
        GR_TMP_INIT_VEC(w, ar, ctx);
        status = gr_csr_mat_mul_vec(w, A, u, ctx);
        _gr_vec_swap(v, w, ar, ctx);
        GR_TMP_CLEAR_VEC(w, ar, ctx);
        return status;
    }

    status = _gr_vec_zero(v, ar, ctx);
    for (i = 0; i < ar; ++i) {
        _gr_csr_mat_borrow_row(row, A, i, ctx);
        status |= gr_sparse_vec_dot_vec(GR_ENTRY(v, i, sz), GR_ENTRY(v, i, sz), 0, row, u, ctx);
    }

    return status;
}

int gr_lil_mat_mul_vec(gr_ptr v, const gr_lil_mat_t A, gr_srcptr u, gr_ctx_t ctx)
{
    slong ar, i, sz;
    int status;

    sz = ctx->sizeof_elem;
    ar = gr_sparse_mat_nrows(A, ctx);

    if (gr_lil_mat_is_zero(A, ctx) == T_TRUE)
    {
        return _gr_vec_zero(v, ar, ctx);
    }

    status = GR_SUCCESS;

    if (u == v)
    {
        gr_ptr w;
        GR_TMP_INIT_VEC(w, ar, ctx);
        status |= gr_lil_mat_mul_vec(w, A, u, ctx);
        _gr_vec_swap(v, w, ar, ctx);
        GR_TMP_CLEAR_VEC(w, ar, ctx);
        return status;
    }

    status |= _gr_vec_zero(v, ar, ctx);
    for (i = 0; i < ar; ++i) {\
        status |= gr_sparse_vec_dot_vec(GR_ENTRY(v, i, sz), GR_ENTRY(v, i, sz), 0, &A->rows[i], u, ctx);
    }

    return status;
}

#define GR_MAT_OOP_FN(status, M, ctx, fn, ...) \
{ \
    gr_mat_t T; \
    gr_mat_init(T, M->r, M->c, ctx); \
    status |= fn(T, __VA_ARGS__, ctx); \
    status |= gr_mat_swap_entrywise(T, M, ctx); \
    gr_mat_clear(T, ctx); \
}

int gr_csr_mat_mul_mat_transpose(gr_mat_t Ct, const gr_csr_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)
{
    slong ar, btr, i;
    int status;

    ar = gr_sparse_mat_nrows(A, ctx);
    btr = gr_mat_nrows(Bt, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_ncols(Bt, ctx) || ar != gr_mat_ncols(Ct, ctx) || btr != gr_mat_nrows(Ct, ctx))
        return GR_DOMAIN;

    if (gr_csr_mat_is_zero(A, ctx) == T_TRUE)
        return gr_mat_zero(Ct, ctx);

    if (Bt == Ct)
    {
        status = GR_SUCCESS;
        GR_MAT_OOP_FN(status, Ct, ctx, gr_csr_mat_mul_mat_transpose, A, Bt);
        return status;
    }

    status = gr_mat_zero(Ct, ctx);
    for (i = 0; i < btr; i++)
    {
        status |= gr_csr_mat_mul_vec(Ct->rows[i], A, Bt->rows[i], ctx);
    }
    return status;
}


int gr_lil_mat_mul_mat_transpose(gr_mat_t Ct, const gr_lil_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)
{
    slong ar, btr, i;
    int status;

    ar = gr_sparse_mat_nrows(A, ctx);
    btr = gr_mat_nrows(Bt, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_ncols(Bt, ctx) || ar != gr_mat_ncols(Ct, ctx) || btr != gr_mat_nrows(Ct, ctx))
        return GR_DOMAIN;

    if (gr_lil_mat_is_zero(A, ctx) == T_TRUE)
        return gr_mat_zero(Ct, ctx);

    if (Bt == Ct)
    {
        status = GR_SUCCESS;
        GR_MAT_OOP_FN(status, Ct, ctx, gr_lil_mat_mul_mat_transpose, A, Bt);
        return status;
    }

    status = gr_mat_zero(Ct, ctx);
    for (i = 0; i < btr; i++) {
        status |= gr_lil_mat_mul_vec(Ct->rows[i], A, Bt->rows[i], ctx);
    }
    return status;
}

int gr_csr_mat_mul_mat(gr_mat_t C, const gr_csr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar, bc, i, j, sz;
    int status;

    ar = gr_sparse_mat_nrows(A, ctx);
    bc = gr_mat_ncols(B, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_nrows(B, ctx) || ar != gr_mat_nrows(C, ctx) || bc != gr_mat_ncols(C, ctx))
        return GR_DOMAIN;

    if (gr_csr_mat_is_zero(A, ctx) == T_TRUE)
        return gr_mat_zero(C, ctx);

    status = GR_SUCCESS;

    if (B == C)
    {
        GR_MAT_OOP_FN(status, C, ctx, gr_csr_mat_mul_mat, A, B);
        return status;
    }

    gr_mat_t Bt;
    gr_mat_t Ct;
    gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);

    sz = ctx->sizeof_elem;

    TMP_INIT;
    TMP_START;
    _GR_MAT_INIT_SHALLOW_TRANSPOSE(Bt, B, sz);
    _GR_MAT_INIT_SHALLOW_TRANSPOSE(Ct, C, sz);

    status |= gr_csr_mat_mul_mat_transpose(Ct, A, Bt, ctx);

    _GR_MAT_SHALLOW_TRANSPOSE(C, Ct, sz);

    flint_free(Bt->rows);
    flint_free(Ct->rows);
    TMP_END;

    return status;
}

int gr_lil_mat_mul_mat(gr_mat_t C, const gr_lil_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar, bc, i, j, sz;
    int status;

    ar = gr_sparse_mat_nrows(A, ctx);
    bc = gr_mat_ncols(B, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_nrows(B, ctx) || ar != gr_mat_nrows(C, ctx) || bc != gr_mat_ncols(C, ctx))
        return GR_DOMAIN;

    if (gr_lil_mat_is_zero(A, ctx) == T_TRUE)
        return gr_mat_zero(C, ctx);

    status = GR_SUCCESS;

    if (B == C)
    {
        GR_MAT_OOP_FN(status, C, ctx, gr_lil_mat_mul_mat, A, B);
        return status;
    }

    gr_mat_t Bt;
    gr_mat_t Ct;
    gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);

    sz = ctx->sizeof_elem;

    TMP_INIT;
    TMP_START;
    _GR_MAT_INIT_SHALLOW_TRANSPOSE(Bt, B, sz);
    _GR_MAT_INIT_SHALLOW_TRANSPOSE(Ct, C, sz);

    status |= gr_lil_mat_mul_mat_transpose(Ct, A, Bt, ctx);

    _GR_MAT_SHALLOW_TRANSPOSE(C, Ct, sz);

    flint_free(Bt->rows);
    flint_free(Ct->rows);
    TMP_END;
    return status;
}
