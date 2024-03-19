/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

void
gr_csr_mat_fit_nnz(gr_csr_mat_t mat, slong nnz, gr_ctx_t ctx)
{
    slong alloc = mat->alloc;
    slong new_alloc = nnz;
    if (new_alloc > alloc)
    {
        slong sz = ctx->sizeof_elem;
        mat->cols = flint_realloc(mat->cols, new_alloc * sizeof(ulong));
        mat->nzs = flint_realloc(mat->nzs, new_alloc * sz);
        _gr_vec_init(GR_ENTRY(mat->nzs, alloc, sz), new_alloc - alloc, ctx);
        mat->alloc = new_alloc;
    }
}

void
gr_coo_mat_fit_nnz(gr_coo_mat_t mat, slong nnz, gr_ctx_t ctx)
{
    slong alloc = mat->alloc;
    slong new_alloc = nnz;
    if (new_alloc > alloc)
    {
        slong sz = ctx->sizeof_elem;
        mat->rows = flint_realloc(mat->rows, new_alloc * sizeof(ulong));
        mat->cols = flint_realloc(mat->cols, new_alloc * sizeof(ulong));
        mat->nzs = flint_realloc(mat->nzs, new_alloc * sz);
        _gr_vec_init(GR_ENTRY(mat->nzs, alloc, sz), new_alloc - alloc, ctx);
        mat->alloc = new_alloc;
    }
}

void
gr_csr_mat_shrink_to_nnz(gr_csr_mat_t mat, gr_ctx_t ctx)
{
    slong nnz = mat->nnz;
    slong sz = ctx->sizeof_elem;
    if (mat->alloc > nnz)
    {
        mat->cols = flint_realloc(mat->cols, nnz * sizeof(ulong));
        _gr_vec_clear(GR_ENTRY(mat->nzs, nnz, sz), mat->alloc - nnz, ctx);
        mat->nzs = flint_realloc(mat->nzs, nnz * sz);
        mat->alloc = nnz;
    }    
}

void
gr_lil_mat_shrink_to_nnz(gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row;

    for (row = 0; row < mat->r; ++row)
    {
        gr_sparse_vec_shrink_to_nnz(&mat->rows[row], ctx);
    }
}

void
gr_coo_mat_shrink_to_nnz(gr_coo_mat_t mat, gr_ctx_t ctx)
{
    slong nnz = mat->nnz;
    slong sz = ctx->sizeof_elem;
    if (mat->alloc > nnz)
    {
        mat->rows = flint_realloc(mat->cols, nnz * sizeof(ulong));
        mat->cols = flint_realloc(mat->cols, nnz * sizeof(ulong));
        _gr_vec_clear(GR_ENTRY(mat->nzs, nnz, sz), mat->alloc - nnz, ctx);
        mat->nzs = flint_realloc(mat->nzs, nnz * sz);
        mat->alloc = nnz;
    }    
}

