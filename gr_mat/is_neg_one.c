/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

truth_t
gr_mat_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    truth_t eq, this_eq;
    slong i, j, r, c, sz;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    if (r == 0 || c == 0)
        return T_TRUE;

    sz = ctx->sizeof_elem;

    eq = T_TRUE;

    for (i = 0; i < r; i++)
    {
        /* todo: use vector functions */
        for (j = 0; j < c; j++)
        {
            if (i == j)
                this_eq = gr_is_neg_one(GR_MAT_ENTRY(mat, i, j, sz), ctx);
            else
                this_eq = gr_is_zero(GR_MAT_ENTRY(mat, i, j, sz), ctx);

            if (this_eq == T_FALSE)
                return T_FALSE;

            if (this_eq == T_UNKNOWN)
                eq = T_UNKNOWN;
        }
    }

    return eq;
}
