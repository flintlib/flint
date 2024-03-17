/*
    Copyright (C) 2024 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_mat_set_csr_mat(gr_mat_t res, const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    ulong row, nz;
    int success = GR_SUCCESS;
    size_t sz = ctx->sizeof_elem;

    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    success |= gr_mat_zero(res, ctx);
    for (row = 0; row < mat->r; ++row)
    {
        for (nz = mat->rows[row]; nz < mat->rows[row + 1]; ++nz)
        {
            success |= gr_set(
                GR_MAT_ENTRY(res, row, mat->cols[nz], sz), 
                GR_ENTRY(mat->entries, nz, sz), 
                ctx
            );
        }
    }
    return success;
}

int gr_mat_set_lil_mat(gr_mat_t res, const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    ulong row;
    int success = GR_SUCCESS;

    if (res->r != mat->r || res->c != mat->c)
        return GR_DOMAIN;

    for (row = 0; row < mat->r; ++row)
    {
        success |= gr_vec_set_sparse_vec(res->rows[row], mat->rows[row], ctx);
    }
    return success;
}
