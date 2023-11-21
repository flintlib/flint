/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int gr_mat_reduce_row(slong * column, gr_mat_t A, slong * P, slong * L, slong m, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = A->c, i, j, r;
    gr_ptr h;
    truth_t is_zero;
    slong sz = ctx->sizeof_elem;

    *column = -1;

    GR_TMP_INIT(h, ctx);

    for (i = 0; i < n; i++)
    {
        is_zero = gr_is_zero(GR_MAT_ENTRY(A, m, i, sz), ctx);

        if (is_zero == T_UNKNOWN)
        {
            status |= GR_UNABLE;
            break;
        }

        if (is_zero == T_FALSE)
        {
            r = P[i];

            if (r != -1)
            {
                for (j = i + 1; j < L[r]; j++)
                {
                    status |= gr_mul(h, GR_MAT_ENTRY(A, r, j, sz), GR_MAT_ENTRY(A, m, i, sz), ctx);
                    status |= gr_sub(GR_MAT_ENTRY(A, m, j, sz), GR_MAT_ENTRY(A, m, j, sz), h, ctx);
                }

                status |= gr_zero(GR_MAT_ENTRY(A, m, i, sz), ctx);
            }
            else
            {
                status |= gr_inv(h, GR_MAT_ENTRY(A, m, i, sz), ctx);
                status |= gr_one(GR_MAT_ENTRY(A, m, i, sz), ctx);

                for (j = i + 1; j < L[m]; j++)
                    status |= gr_mul(GR_MAT_ENTRY(A, m, j, sz), GR_MAT_ENTRY(A, m, j, sz), h, ctx);

                P[i] = m;

                *column = i;
                break;
            }
        }
    }

    GR_TMP_CLEAR(h, ctx);

    return status;
}
