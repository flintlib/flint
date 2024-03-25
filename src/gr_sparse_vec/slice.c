/*
    Copyright (C) 2024 Kartik Venkatram and Alden Walker

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_sparse_vec.h"

int
gr_sparse_vec_slice(gr_sparse_vec_t dst, const gr_sparse_vec_t src, slong col_start, slong col_end, gr_ctx_t ctx)
{
    slong sz,i,nnz, new_nnz, i_start, i_end, status;
    nnz = src->nnz;
    sz = ctx->sizeof_elem;
    i_start = 0;
    i_end = nnz;
    status = GR_SUCCESS;
    for (i = 0; i < nnz; i++)
    {
        /* If we find a valid column, start the interval */
        if (src->inds[i] >= col_start)
        {
            i_start = i;
            break;
        }
    }
    for (; i <= nnz; i++)
    {
        /* Note we will always do this; the first time we see a column outside the interval, or we run off the end, stop it */
        if (i == nnz || src->inds[i] >= col_end)
        {
            i_end = i;
            break;
        }
    }
    /* We are messing with the internals of dst; this should be ok even if it is also src */
    new_nnz = i_end - i_start;
    dst->length = col_end - col_start;
    gr_sparse_vec_fit_nnz(dst, new_nnz, ctx);
    for (i = i_start; i < i_end; i++)
    {
        dst->inds[i-i_start] = dst->inds[i] - col_start;
        status |= gr_set(GR_ENTRY(dst->nzs, i-i_start, sz), GR_ENTRY(src->nzs, i, sz), ctx);
    }
    dst->nnz = new_nnz;
    return status;
}
