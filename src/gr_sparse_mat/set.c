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
    gr_sparse_vec_t res_row;

    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    gr_csr_mat_fit_nnz(res, mat->nnz, ctx);

    res->rows[0] = 0;
    for(row = 0; row < mat->r; row++) {
        res->rows[row+1] = res->rows[row] + mat->rows[row]->nnz;
        _gr_csr_mat_borrow_row(res_row, res, row, ctx);
        success |= gr_sparse_vec_set(res_row, mat->rows[row], ctx);
    }
    res->nnz = mat->nnz;

    return success;
}

int gr_lil_mat_set_csr_mat(gr_lil_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    ulong row;
    int success = GR_SUCCESS;
    gr_sparse_vec_t mat_row;
    
    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    for(row = 0; row < mat->r; row++) {
        _gr_csr_mat_borrow_row(mat_row, mat, row, ctx);
        success |= gr_sparse_vec_set(res->rows[row], mat_row, ctx);
    }
    res->nnz = mat->nnz;

    return success;
}
