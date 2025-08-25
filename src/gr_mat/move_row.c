/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "gr.h"
#include "gr_mat.h"

int gr_mat_move_row(gr_mat_t A, slong i, slong new_i, gr_ctx_t ctx)
{
    gr_ptr tmp;
    slong c = A->c, r = A->c;
    slong sz = ctx->sizeof_elem;
    slong j;

    if (i < 0 || i >= r || new_i < 0 || i >= r)
        return GR_DOMAIN;

    if (i == new_i)
        return GR_SUCCESS;

    tmp = GR_TMP_ALLOC(c * sz);

    if (new_i < i)
    {
        memcpy(tmp, GR_MAT_ENTRY(A, i, 0, sz), c * sz);
        for (j = i; j > new_i; j--)
            memcpy(GR_MAT_ENTRY(A, j, 0, sz), GR_MAT_ENTRY(A, j - 1, 0, sz), c * sz);
        memcpy(GR_MAT_ENTRY(A, new_i, 0, sz), tmp, c * sz);
    }
    else
    {
        memcpy(tmp, GR_MAT_ENTRY(A, i, 0, sz), c * sz);
        for (j = i; j < new_i; j++)
            memcpy(GR_MAT_ENTRY(A, j, 0, sz), GR_MAT_ENTRY(A, j + 1, 0, sz), c * sz);
        memcpy(GR_MAT_ENTRY(A, new_i, 0, sz), tmp, c * sz);
    }

    GR_TMP_FREE(tmp, c * sz);

    return GR_SUCCESS;
}
