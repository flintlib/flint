/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_mat.h"

truth_t gr_csr_mat_is_one(const gr_csr_mat_t mat, gr_ctx_t ctx) 
{
    slong row, idx;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    gr_method_unary_predicate is_one = GR_UNARY_PREDICATE(ctx, IS_ONE);
    truth_t this_eq;
    truth_t ret = T_TRUE;
    slong sz = ctx->sizeof_elem;
    
    for (row = 0; row < mat->r; ++row)
    {
        for (idx = mat->rows[row]; idx < mat->rows[row+1]; idx++)
        {
            this_eq = (mat->cols[idx] == row ? is_one : is_zero)(GR_ENTRY(mat->entries, idx, sz), ctx);
            if (this_eq == T_FALSE)
                return T_FALSE;
            else
                ret = T_UNKNOWN;
        }
    }
    return ret;
}

truth_t gr_lil_mat_is_one(const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    slong row, idx;
    gr_method_unary_predicate is_zero = GR_UNARY_PREDICATE(ctx, IS_ZERO);
    gr_method_unary_predicate is_one = GR_UNARY_PREDICATE(ctx, IS_ONE);
    truth_t this_eq;
    truth_t ret = T_TRUE;
    slong sz = ctx->sizeof_elem;
    
    for (row = 0; row < mat->r; ++row)
    {
        for (idx = 0; idx < mat->rows[row]->nnz; idx++)
        {
            this_eq = (mat->rows[row]->inds[idx] == row ? is_one : is_zero)(GR_ENTRY(mat->rows[row]->entries, idx, sz), ctx);
            if (this_eq == T_FALSE)
                return T_FALSE;
            else
                ret = T_UNKNOWN;
        }
    }
    return ret;
}
