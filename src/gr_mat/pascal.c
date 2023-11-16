/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: optimize for fmpz */

int
gr_mat_pascal(gr_mat_t mat, int triangular, gr_ctx_t ctx)
{
    slong R, C, i, j;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

#define ENTRY(i, j) GR_MAT_ENTRY(mat, i, j, sz)

    R = gr_mat_nrows(mat, ctx);
    C = gr_mat_ncols(mat, ctx);

    if (R == 0 || C == 0)
        return status;

    if (triangular == 1)
    {
        for (i = 1; i < R; i++)
            for (j = 0; j < i && j < C; j++)
                status |= gr_zero(ENTRY(i, j), ctx);

        for (j = 0; j < C; j++)
            status |= gr_one(ENTRY(0, j), ctx);

        for (i = 1; i < R && i < C; i++)
            status |= gr_one(ENTRY(i, i), ctx);

        for (i = 1; i < R; i++)
            for (j = i + 1; j < C; j++)
                status |= gr_add(ENTRY(i, j),
                    ENTRY(i, j - 1), ENTRY(i - 1, j - 1), ctx);
    }
    else if (triangular == -1)
    {
        for (i = 0; i < R; i++)
            for (j = i + 1; j < C; j++)
                status |= gr_zero(ENTRY(i, j), ctx);

        for (i = 0; i < R; i++)
            status |= gr_one(ENTRY(i, 0), ctx);

        for (i = 1; i < R && i < C; i++)
            status |= gr_one(ENTRY(i, i), ctx);

        for (i = 2; i < R; i++)
            for (j = 1; j < i && j < C; j++)
                status |= gr_add(ENTRY(i, j),
                    ENTRY(i - 1, j - 1), ENTRY(i - 1, j), ctx);
    }
    else if (triangular == 0)
    {
        for (j = 0; j < C; j++)
            status |= gr_one(ENTRY(0, j), ctx);

        for (i = 1; i < R; i++)
            status |= gr_one(ENTRY(i, 0), ctx);

        for (i = 1; i < R; i++)
            for (j = 1; j < C; j++)
                status |= gr_add(ENTRY(i, j),
                    ENTRY(i, j - 1), ENTRY(i - 1, j), ctx);
    }
    else
    {
        status = GR_DOMAIN;
    }

    return status;
}
