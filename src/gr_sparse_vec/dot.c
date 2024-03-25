/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "gr_sparse_vec.h"

int gr_sparse_vec_dot(gr_ptr dst, gr_srcptr initial, int subtract, const gr_sparse_vec_t src1, const gr_sparse_vec_t src2, gr_ctx_t ctx)
{
    int status;
    slong nz_idx1, nz_idx2;
    slong sz = ctx->sizeof_elem;
    
    if (src1->length != src2->length)
    {
        return GR_DOMAIN;
    }
    status = gr_set(dst, initial, ctx);
    for (nz_idx1 = 0, nz_idx2 = 0; nz_idx1 < src1->nnz && nz_idx2 < src2->nnz; )
    {
        if (src1->inds[nz_idx1] < src2->inds[nz_idx2])
            nz_idx1++;
        else if (src1->inds[nz_idx1] > src2->inds[nz_idx2])
            nz_idx2++;
        else {
            if (subtract)
                status |= gr_submul(dst, GR_ENTRY(src1->nzs, nz_idx1, sz), GR_ENTRY(src2->nzs, nz_idx2, sz), ctx);
            else
                status |= gr_addmul(dst, GR_ENTRY(src1->nzs, nz_idx1, sz), GR_ENTRY(src2->nzs, nz_idx2, sz), ctx);
            nz_idx1++, nz_idx2++;
        }
    }
    return status;
}

int gr_sparse_vec_dot_vec(gr_ptr dst, gr_srcptr initial, int subtract, const gr_sparse_vec_t src1, gr_srcptr src2, gr_ctx_t ctx)
{
    int status;
    slong nz_idx;
    slong sz = ctx->sizeof_elem;
    
    status = gr_set(dst, initial, ctx);
    //flint_printf("dst = "); status |= gr_println(dst, ctx);
    for (nz_idx = 0; nz_idx < src1->nnz; nz_idx++)
    {
        if (subtract)
            status |= gr_submul(dst, GR_ENTRY(src1->nzs, nz_idx, sz), GR_ENTRY(src2, src1->inds[nz_idx], sz), ctx);
        else
            status |= gr_addmul(dst, GR_ENTRY(src1->nzs, nz_idx, sz), GR_ENTRY(src2, src1->inds[nz_idx], sz), ctx);
        //flint_printf("dst = "); status |= gr_println(dst, ctx);
    }
    return status;
}
