/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_mat.h"

int
gr_mat_transpose(gr_mat_t B, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (B->r != A->c || B->c != A->r)
        return GR_DOMAIN;

    if (A->r == 0 || A->c == 0)
        return GR_SUCCESS;

    if (A == B)  /* In-place, guaranteed to be square */
    {
        slong stride = A->stride;

        /* Optimize for common sizes */
        if (ctx->sizeof_elem == sizeof(ulong))
        {
            ulong * a = A->entries;

            for (i = 0; i < A->r - 1; i++)
                for (j = i + 1; j < A->c; j++)
                    FLINT_SWAP(ulong, a[i * stride + j], a[j * stride + i]);
        }
        else
        {
            gr_ptr a = A->entries;

            for (i = 0; i < A->r - 1; i++)
                for (j = i + 1; j < A->c; j++)
                    gr_swap(GR_ENTRY(a, i * stride + j, sz), GR_ENTRY(a, j * stride + i, sz), ctx);
        }
    }
    else  /* Not aliased; general case */
    {
        slong Astride = A->stride;
        slong Bstride = B->stride;
        gr_srcptr a = A->entries;
        gr_ptr b = B->entries;
        gr_method_unary_op set = GR_UNARY_OP(ctx, SET);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                status |= set(GR_ENTRY(b, j * Bstride + i, sz), GR_ENTRY(a, i * Astride + j, sz), ctx);
    }

    return status;
}
