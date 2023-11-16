/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

truth_t
gr_mat_is_scalar(const gr_mat_t mat, gr_ctx_t ctx)
{
    truth_t eq, this_eq;
    slong i, ar, ac, sz;

    ar = gr_mat_nrows(mat, ctx);
    ac = gr_mat_ncols(mat, ctx);

    sz = ctx->sizeof_elem;

    eq = gr_mat_is_diagonal(mat, ctx);

    if (eq != T_FALSE)
    {
        for (i = 1; i < FLINT_MIN(ar, ac); i++)
        {
            this_eq = gr_equal(GR_MAT_ENTRY(mat, i, i, sz), GR_MAT_ENTRY(mat, 0, 0, sz), ctx);

            if (this_eq == T_FALSE)
                return T_FALSE;

            if (this_eq == T_UNKNOWN)
                eq = T_UNKNOWN;
        }
    }

    return eq;
}
