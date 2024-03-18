/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

truth_t
gr_csr_mat_equal(const gr_csr_mat_t mat1, const gr_csr_mat_t mat2, gr_ctx_t ctx)
{
    slong i, r, c, nnz;

    r = gr_sparse_mat_nrows(mat1, ctx);
    c = gr_sparse_mat_ncols(mat1, ctx);
    nnz = gr_sparse_mat_nnz(mat1, ctx);

    if (r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx) ||
        nnz != gr_sparse_mat_nnz(mat2, ctx))
    {
        return T_FALSE;
    }

    if (r == 0 || c == 0)
        return T_TRUE;

    // TODO: can use memcmp?
    for (i = 0; i < nnz; ++i)
    {
        if (mat1->cols[i] != mat2->cols[i])
            return T_FALSE;
    }

    return _gr_vec_equal(
        mat1->entries,
        mat2->entries,
        nnz,
        ctx
    );
}

truth_t
gr_lil_mat_equal(const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{
    slong i, r, c, nnz;
    truth_t row_is_eq;

    r = gr_sparse_mat_nrows(mat1, ctx);
    c = gr_sparse_mat_ncols(mat1, ctx);
    nnz = gr_sparse_mat_nnz(mat1, ctx);

    if (r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx) ||
        nnz != gr_sparse_mat_nnz(mat2, ctx))
    {
        return T_FALSE;
    }

    if (r == 0 || c == 0)
        return T_TRUE;

    for (i = 0; i < r; ++i)
    {
        row_is_eq = gr_sparse_vec_equal(mat1->rows[i], mat2->rows[i], ctx);
        if (row_is_eq != T_TRUE)
            break;
    }

    return row_is_eq;
}
