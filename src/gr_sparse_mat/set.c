/*
    Copyright (C) 2024 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_csr_mat_set(gr_csr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx) {
    int status;

    if (res->r != mat->r || res->c != mat->c) {
        return GR_DOMAIN;
    }
    res->nnz = mat->nnz;
    gr_csr_mat_fit_nnz(res, mat->nnz, ctx);
    memcpy(res->rows, mat->rows, mat->r * sizeof(ulong));
    memcpy(res->cols, mat->cols, mat->nnz * sizeof(ulong));
    status = _gr_vec_set(res->entries, mat->entries, mat->nnz, ctx);
    return status;
}

int gr_lil_mat_set(gr_lil_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx) {
    ulong row;
    int success = GR_SUCCESS;

    if (res->r != mat->r || res->c != mat->c) {
        return GR_DOMAIN;
    }
    res->nnz = mat->nnz;

    for (row = 0; row < mat->r; ++row) {
        success |= gr_sparse_vec_set(res->rows[row], mat->rows[row], ctx);
    }
    return success;
}

int gr_csr_mat_set_lil_mat(gr_csr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    ulong row;
    int success = GR_SUCCESS;
    gr_sparse_vec_struct *mat_row;

    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    gr_csr_mat_fit_nnz(res, mat->nnz, ctx);

    res->nnz = 0;
    for(row = 0; row < mat->r; row++) {
        mat_row = mat->rows[row];
        res->rows[row] = res->nnz;
        memcpy(res->cols + res->nnz, mat_row->cols, mat_row->nnz);
        success |= _gr_vec_set(GR_ENTRY(res->entries, res->nnz, ctx->sizeof_elem), mat_row->entries, mat_row->nnz, ctx);
        res->nnz += mat_row->nnz;
    }
    res->rows[res->r] = res->nnz;

    return success;
}

int gr_lil_mat_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    ulong row;
    int success = GR_SUCCESS;
    gr_sparse_vec_t tmp_vec;
    
    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    tmp_vec->alloc = 0; // All memory borrowed
    tmp_vec->length = mat->c;
    for(row = 0; row < mat->r; row++) {
        tmp_vec->cols = mat->cols + mat->rows[row];
        tmp_vec->entries = GR_ENTRY(mat->entries, mat->rows[row], ctx->sizeof_elem);
        tmp_vec->nnz = mat->rows[row+1] - mat->rows[row];
        success |= gr_sparse_vec_set(res->rows[row], tmp_vec, ctx);
    }
    res->nnz = mat->nnz;

    return success;
}
