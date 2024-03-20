/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gr_sparse_mat.h"

int gr_csr_mat_permute_cols(gr_csr_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    slong row;
    gr_sparse_vec_t tmp;
    int status = GR_SUCCESS;

    for (row = 0; row < mat->r; ++row)
    {
        _gr_csr_mat_borrow_row(tmp, mat, row, ctx);
        status |= gr_sparse_vec_permute_inds(tmp, tmp, perm, ctx);
    }
    return status;
}

int gr_lil_mat_permute_cols(gr_lil_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    slong row;
    int status = GR_SUCCESS;

    for (row = 0; row < mat->r; ++row)
    {
        status |= gr_sparse_vec_permute_inds(&mat->rows[row], &mat->rows[row], perm, ctx);
    }
    return status;
}

int gr_coo_mat_permute_cols(gr_coo_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    slong nz_idx;
    int status = GR_SUCCESS;

    for (nz_idx = 0; nz_idx < mat->nnz; ++nz_idx)
    {
        mat->cols[nz_idx] = perm[mat->cols[nz_idx]];
    }
    if (mat->is_canonical)
    {
        gr_coo_mat_canonicalize(mat, ctx);
    }
    return status;
}