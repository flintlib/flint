/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_mul_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar, ac, br, bc, i, j, sz;
    int status;

    ar = gr_mat_nrows(A, ctx);
    ac = gr_mat_ncols(A, ctx);
    br = gr_mat_nrows(B, ctx);
    bc = gr_mat_ncols(B, ctx);

    if (ac != br || ar != gr_mat_nrows(C, ctx) || bc != gr_mat_ncols(C, ctx))
        return GR_DOMAIN;

    if (br == 0)
    {
        return gr_mat_zero(C, ctx);
    }

    status = GR_SUCCESS;

    if (A == C || B == C)
    {
        gr_mat_t T;
        gr_mat_init(T, ar, bc, ctx);
        status |= gr_mat_mul_classical(T, A, B, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    sz = ctx->sizeof_elem;

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                status |= gr_mul(GR_MAT_ENTRY(C, i, j, sz),
                                 GR_MAT_ENTRY(A, i, 0, sz),
                                 GR_MAT_ENTRY(B, 0, j, sz), ctx);
            }
        }
    }
    else
    {
        gr_ptr tmp;
        gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sz * br * bc);

        /* Make a shallow transpose so that we can use dot products.
           Inline common sizes. (Caution: are we sure about the alignment?
           Some asserts would be nice here.)
           Todo: we may want inlining in nonsingular_solve etc. as well. */
        for (i = 0; i < br; i++)
        {
            for (j = 0; j < bc; j++)
            {
                switch (sz)
                {
#if 0
                    case 1:
                        ((int8_t *) GR_ENTRY(tmp, j * br + i, 1))[0] = ((int8_t *) GR_MAT_ENTRY(B, i, j, 1))[0];
                        break;
                    case 2:
                        ((int16_t *) GR_ENTRY(tmp, j * br + i, 2))[0] = ((int16_t *) GR_MAT_ENTRY(B, i, j, 2))[0];
                        break;
                    case 4:
                        ((int32_t *) GR_ENTRY(tmp, j * br + i, 4))[0] = ((int32_t *) GR_MAT_ENTRY(B, i, j, 4))[0];
                        break;
#if FLINT_BITS == 64
                    case 8:
                        ((int64_t *) GR_ENTRY(tmp, j * br + i, 8))[0] = ((int64_t *) GR_MAT_ENTRY(B, i, j, 8))[0];
                        break;
#endif
#endif
                    default:
                        set_shallow(GR_ENTRY(tmp, j * br + i, sz), GR_MAT_ENTRY(B, i, j, sz), ctx);
                }
            }
        }

        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                status |= _gr_vec_dot(GR_MAT_ENTRY(C, i, j, sz), NULL, 0,
                    GR_MAT_ENTRY(A, i, 0, sz), GR_ENTRY(tmp, j * br, sz), br, ctx);
            }
        }

        TMP_END;
    }

    return status;
}
