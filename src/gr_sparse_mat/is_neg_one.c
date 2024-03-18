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
gr_csr_mat_is_neg_one(const gr_csr_mat_t mat, gr_ctx_t ctx)
{
    int row;
    gr_method_unary_predicate is_neg_one = GR_UNARY_PREDICATE(ctx, IS_NEG_ONE);
    truth_t this_eq;
    slong sz = ctx->sizeof_elem;
    
    if (mat->r != mat->nnz)
        return T_FALSE;

    for (row = 0; row < mat->r; ++row)
    {
        if (mat->cols[row] != row)
            return T_FALSE;

        this_eq = is_neg_one(GR_ENTRY(mat->entries, row, sz), ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;
        else if (this_eq == T_UNKNOWN)
            return T_UNKNOWN;
    }
    return T_TRUE;
}

truth_t gr_lil_mat_is_neg_one(const gr_lil_mat_t mat, gr_ctx_t ctx)
{
    int row;
    gr_method_unary_predicate is_neg_one = GR_UNARY_PREDICATE(ctx, IS_NEG_ONE);
    truth_t this_eq;
    
    if (mat->r != mat->nnz)
        return T_FALSE;

    for (row = 0; row < mat->r; ++row)
    {
        if (mat->rows[row]->nnz != 1 || mat->rows[row]->inds[0] != row)
            return T_FALSE;

        this_eq = is_neg_one(mat->rows[row]->entries, ctx);

        if (this_eq == T_FALSE)
            return T_FALSE;
        else if (this_eq == T_UNKNOWN)
            return T_UNKNOWN;
    }
    return T_TRUE;
}
