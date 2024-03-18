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
    slong row;
    gr_sparse_vec_t tmp1, tmp2;
    truth_t row_is_eq;
    truth_t ret = T_TRUE;

    if (gr_mat_is_compatible(mat1, mat2, ctx) == T_FALSE)
    {
        return T_FALSE;
    }
    for (row = 0; row < mat1->r; row++)
    {
        _gr_csr_mat_borrow_row(tmp1, mat1, row, ctx);
        _gr_csr_mat_borrow_row(tmp2, mat2, row, ctx);
        row_is_eq = gr_sparse_vec_equal(tmp1, tmp2, ctx);
        if (row_is_eq == T_FALSE)
            return T_FALSE;
        else if (row_is_eq == T_UNKNOWN)
            ret = T_UNKNOWN;
    }
    return ret;
}

truth_t
gr_lil_mat_equal(const gr_lil_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{
    slong row;
    truth_t row_is_eq;
    truth_t ret = T_TRUE;

    if (gr_mat_is_compatible(mat1, mat2, ctx) == T_FALSE)
    {
        return T_FALSE;
    }
    for (row = 0; row < mat1->r; row++)
    {
        row_is_eq = gr_sparse_vec_equal(mat1->rows[row], mat2->rows[row], ctx);
        if (row_is_eq == T_FALSE)
            return T_FALSE;
        else if (row_is_eq == T_UNKNOWN)
            ret = T_UNKNOWN;
    }
    return ret;
}

truth_t
gr_csr_mat_equal_lil_mat(const gr_csr_mat_t mat1, const gr_lil_mat_t mat2, gr_ctx_t ctx)
{
    slong row;
    gr_sparse_vec_t tmp;
    truth_t row_is_eq;
    truth_t ret = T_TRUE;

    if (gr_mat_is_compatible(mat1, mat2, ctx) == T_FALSE)
    {
        return T_FALSE;
    }
    for (row = 0; row < mat1->r; row++)
    {
        _gr_csr_mat_borrow_row(tmp, mat1, row, ctx);
        row_is_eq = gr_sparse_vec_equal(tmp, mat2->rows[row], ctx);
        if (row_is_eq == T_FALSE)
            return T_FALSE;
        else if (row_is_eq == T_UNKNOWN)
            ret = T_UNKNOWN;
    }
    return ret;
}

