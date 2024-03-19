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

int gr_csr_mat_mul_vec(gr_vec_t v, const gr_csr_mat_t A, const gr_vec_t u, gr_ctx_t ctx)
{
    slong ar, i;
    int status;
    gr_sparse_vec_t row;
    gr_ptr res;

    ar = gr_sparse_mat_nrows(A, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_vec_length(u, ctx) || ar != gr_vec_length(v, ctx))
        return GR_DOMAIN;

    if (gr_csr_mat_is_zero(A, ctx) == T_TRUE)
    {
        return _gr_vec_zero(v->entries, v->length, ctx);
    }

    status = GR_SUCCESS;

    if (u == v)
    {
        gr_vec_t w;
        gr_vec_init(w, ar, ctx);
        status |= gr_csr_mat_mul_vec(w, A, u, ctx);
        _gr_vec_swap(v, w, ar, ctx);
        gr_vec_clear(w, ctx);
        return status;
    }

    status |= _gr_vec_zero(v->entries, v->length, ctx);
    for (i = 0; i < ar; ++i) {
        _gr_csr_mat_borrow_row(row, A, i, ctx);
        res = gr_vec_entry_ptr(v, i, ctx);
        status |= gr_sparse_vec_dot_vec(res, res, 0, row, u, ctx);
    }

    return status;
}

int gr_lil_mat_mul_vec(gr_vec_t v, const gr_lil_mat_t A, const gr_vec_t u, gr_ctx_t ctx)
{
    slong ar, i;
    int status;
    gr_ptr res;

    ar = gr_sparse_mat_nrows(A, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_vec_length(u, ctx) || ar != gr_vec_length(v, ctx))
        return GR_DOMAIN;


    if (gr_lil_mat_is_zero(A, ctx) == T_TRUE)
    {
        return _gr_vec_zero(v->entries, v->length, ctx);
    }

    status = GR_SUCCESS;

    if (u == v)
    {
        gr_vec_t w;
        gr_vec_init(w, ar, ctx);
        status |= gr_lil_mat_mul_vec(w, A, u, ctx);
        _gr_vec_swap(v, w, ar, ctx);
        gr_vec_clear(w, ctx);
        return status;
    }

    status |= _gr_vec_zero(v->entries, v->length, ctx);
    for (i = 0; i < ar; ++i) {
        res = gr_vec_entry_ptr(v, i, ctx);
        status |= gr_sparse_vec_dot_vec(res, res, 0, &A->rows[i], u, ctx);
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

void _gr_mat_borrow_row(gr_vec_t row, const gr_mat_t mat, slong i, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    row->length = mat->c;
    row->entries = GR_MAT_ENTRY(mat, i, 0, sz);
}


int gr_csr_mat_mul_mat_transpose(gr_mat_t Ct, const gr_csr_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)
{
    slong ar, btr, i;
    int status;
    gr_vec_t bt_row, ct_row;

    ar = gr_sparse_mat_nrows(A, ctx);
    btr = gr_mat_nrows(Bt, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_ncols(Bt, ctx) || ar != gr_mat_ncols(Ct, ctx) || btr != gr_mat_nrows(Ct, ctx))
        return GR_DOMAIN;

    status = gr_mat_zero(Ct, ctx);

    if (gr_csr_mat_is_zero(A, ctx))
    {
        return status;
    }

    if (Bt == Ct)
    {
        GR_MAT_OOP_FN(status, Ct, ctx, gr_csr_mat_mul_mat_transpose, A, Bt);
        return status;
    }

    for (i = 0; i < btr; i++)
    {
        _gr_mat_borrow_row(ct_row, Ct, i, ctx);
        _gr_mat_borrow_row(bt_row, Bt, i, ctx);
        status |= gr_csr_mat_mul_vec(ct_row, A, bt_row, ctx);
    }
    return status;
}


int gr_lil_mat_mul_mat_transpose(gr_mat_t Ct, const gr_lil_mat_t A, const gr_mat_t Bt, gr_ctx_t ctx)
{
    slong ar, btr, i;
    int status;
    gr_vec_t bt_row, ct_row;

    ar = gr_sparse_mat_nrows(A, ctx);
    btr = gr_mat_nrows(Bt, ctx);

    if (gr_sparse_mat_ncols(A, ctx) != gr_mat_ncols(Bt, ctx) || ar != gr_mat_ncols(Ct, ctx) || btr != gr_mat_nrows(Ct, ctx))
        return GR_DOMAIN;

    status = gr_mat_zero(Ct, ctx);

    if (gr_lil_mat_is_zero(A, ctx))
    {
        return status;
    }

    status = GR_SUCCESS;

    if (Bt == Ct)
    {
        GR_MAT_OOP_FN(status, Ct, ctx, gr_lil_mat_mul_mat_transpose, A, Bt);
        return status;
    }

    for (i = 0; i < btr; i++) {
        _gr_mat_borrow_row(ct_row, Ct, i, ctx);
        _gr_mat_borrow_row(bt_row, Bt, i, ctx);
        status |= gr_lil_mat_mul_vec(ct_row, A, bt_row, ctx);
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

    if (gr_csr_mat_is_zero(A, ctx))
    {
        return gr_mat_zero(C, ctx);
    }

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

    if (gr_lil_mat_is_zero(A, ctx))
    {
        return gr_mat_zero(C, ctx);
    }

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

    flint_free(Bt->rows);
    flint_free(Ct->rows);
    TMP_END;
    return status;
}
