/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "gr_mat.h"

int
gr_mat_set_fmpz_mat(gr_mat_t res, const fmpz_mat_t mat, gr_ctx_t ctx)
{
    slong i, j;
    slong m = mat->r;
    slong n = mat->c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            status |= gr_set_fmpz(GR_MAT_ENTRY(res, i, j, sz), fmpz_mat_entry(mat, i, j), ctx);
        }
    }

    return status;
}
