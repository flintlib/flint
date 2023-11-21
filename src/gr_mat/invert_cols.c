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
gr_mat_invert_cols(gr_mat_t mat, slong * perm, gr_ctx_t ctx)
{
    if (gr_mat_is_empty(mat, ctx) == T_FALSE)
    {
        slong sz = ctx->sizeof_elem;
        slong t, i;
        slong c = mat->c;
        slong k = mat->c / 2;

        if (perm != NULL)
            for (i = 0; i < k; i++)
                FLINT_SWAP(slong, perm[i], perm[c - i - 1]);

        for (t = 0; t < mat->r; t++)
            for (i = 0; i < k; i++)
                gr_swap(GR_MAT_ENTRY(mat, t, i, sz), GR_MAT_ENTRY(mat, t, c - i - 1, sz), ctx);
    }

    return GR_SUCCESS;
}
