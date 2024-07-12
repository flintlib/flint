/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

int gr_csr_mat_invert_rows(gr_csr_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong r = mat->r;
    slong nnz = mat->nnz;
    slong i, j, k;
    slong sz = ctx->sizeof_elem;

    // Handle permutation if provided
    if (perm != NULL)
    {
        for (i = 0; i < r / 2; i++)
        {
            FLINT_SWAP(slong, perm[i], perm[r - i - 1]);
        }
    }

    // Reverse row offsets
    for (i = 0; i < r / 2; ++i)
    {
        FLINT_SWAP(ulong, mat->rows[i + 1], mat->rows[r - i]);
    }

    // Reverse all columns and elements
    for (j = 0; j < nnz / 2; ++j)
    {
        k = nnz - j - 1;
        FLINT_SWAP(ulong, mat->cols[j], mat->cols[k]);
        gr_swap(GR_ENTRY(mat->nzs, j, sz), GR_ENTRY(mat->nzs, k, sz), ctx);
    }

    // Reverse columns and elements in each row
    for (i = 0; i < r; ++i)
    {
        // Fix row offset
        mat->rows[i+1] -= (i < r - 1 ? mat->rows[i+2] : 0) - mat->rows[i];
        for (j = mat->rows[i]; j < mat->rows[i+1] / 2; ++j)
        {
            k = mat->rows[i+1] - j - 1;
            FLINT_SWAP(ulong, mat->cols[j], mat->cols[k]);
            gr_swap(GR_ENTRY(mat->nzs, j, sz), GR_ENTRY(mat->nzs, k, sz), ctx);
        }
    }
    return status;
}

int gr_lil_mat_invert_rows(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong r = mat->r;
    slong i;

    for (i = 0; i < r / 2; i++)
        status |= gr_lil_mat_swap_rows(mat, perm, i, r - i - 1, ctx);

    return status;
}
