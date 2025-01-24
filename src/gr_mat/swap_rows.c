/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int gr_mat_swap_rows(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
{
    /* todo: bounds checking */

    if (r != s && gr_mat_is_empty(mat, ctx) == T_FALSE)
    {
        slong sz = ctx->sizeof_elem;

        if (perm != NULL)
            FLINT_SWAP(slong, perm[r], perm[s]);

        _gr_vec_swap(GR_MAT_ENTRY(mat, r, 0, sz), GR_MAT_ENTRY(mat, s, 0, sz), mat->c, ctx);
    }

    return GR_SUCCESS;
}
