/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_swap_cols(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
{
    /* todo: bounds checking */

    if (r != s && gr_mat_is_empty(mat, ctx) == T_FALSE)
    {
        slong t;
        slong sz = ctx->sizeof_elem;

        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        for (t = 0; t < mat->r; t++)
            gr_swap(GR_MAT_ENTRY(mat, t, r, sz), GR_MAT_ENTRY(mat, t, s, sz), ctx);
    }

    return GR_SUCCESS;
}
