/*
    Copyright (C) 2024 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_mat_set_csr_mat(gr_mat_t dst, const gr_csr_mat_t src, gr_ctx_t ctx)
{
    ulong row, nz;
    int status = GR_SUCCESS;
    size_t sz = ctx->sizeof_elem;

    if (dst->r != src->r || dst->c != src->c)
        return GR_DOMAIN;

    status |= gr_mat_zero(dst, ctx);
    for (row = 0; row < src->r; ++row)
    {
        for (nz = src->rows[row]; nz < src->rows[row + 1]; ++nz)
        {
            status |= gr_set(
                GR_MAT_ENTRY(dst, row, src->cols[nz], sz), 
                GR_ENTRY(src->nzs, nz, sz), 
                ctx
            );
        }
    }
    return status;
}

int gr_mat_set_lil_mat(gr_mat_t dst, const gr_lil_mat_t src, gr_ctx_t ctx)
{
    ulong row;
    int status = GR_SUCCESS;

    if (dst->r != src->r || dst->c != src->c)
        return GR_DOMAIN;

    for (row = 0; row < src->r; ++row)
    {
        status |= gr_vec_set_sparse_vec(dst->rows[row], &src->rows[row], ctx);
    }
    return status;
}

int gr_mat_set_coo_mat(gr_mat_t dst, const gr_coo_mat_t src, gr_ctx_t ctx)
{
    ulong nz;
    int status = GR_SUCCESS;
    gr_ptr dst_entry;
    gr_srcptr src_entry;
    size_t sz = ctx->sizeof_elem;

    if (dst->r != src->r || dst->c != src->c)
        return GR_DOMAIN;

    status |= gr_mat_zero(dst, ctx);
    for (nz = 0; nz < src->nnz; ++nz)
    {
        dst_entry = GR_MAT_ENTRY(dst, src->rows[nz], src->cols[nz], sz);
        src_entry = GR_ENTRY(src->nzs, nz, sz);
        if (src->is_canonical)
            status |= gr_set(dst_entry, src_entry, ctx);
        else
            status |= gr_add(dst_entry, dst_entry, src_entry, ctx);
    }
    return status;
}
