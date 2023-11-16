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
gr_mat_find_nonzero_pivot_large_abs(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    slong best_row, i;
    int unknown;
    int cmp, comp_ok;
    truth_t is_zero;
    slong sz;

    if (end_row <= start_row)
        return GR_DOMAIN;

    best_row = -1;
    unknown = 0;

    sz = ctx->sizeof_elem;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = gr_is_zero(GR_MAT_ENTRY(mat, i, column, sz), ctx);

        if (is_zero == T_FALSE)
        {
            if (best_row == -1)
            {
                best_row = i;
            }
            else
            {
                comp_ok = gr_cmpabs(&cmp, GR_MAT_ENTRY(mat, i, column, sz), GR_MAT_ENTRY(mat, best_row, column, sz), ctx);

                if (comp_ok == GR_SUCCESS && cmp > 0)
                    best_row = i;
            }
        }

        if (is_zero == T_UNKNOWN)
            unknown = 1;
    }

    if (best_row == -1)
    {
        *pivot_row = -1;
        if (unknown)
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }
    else
    {
        *pivot_row = best_row;
        return GR_SUCCESS;
    }
}

int
gr_mat_find_nonzero_pivot_generic(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    slong i;
    int unknown;
    truth_t is_zero;
    slong sz;

    if (end_row <= start_row)
        return GR_DOMAIN;

    unknown = 0;
    sz = ctx->sizeof_elem;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = gr_is_zero(GR_MAT_ENTRY(mat, i, column, sz), ctx);

        if (is_zero == T_FALSE)
        {
            *pivot_row = i;
            return GR_SUCCESS;
        }

        if (is_zero == T_UNKNOWN)
            unknown = 1;
    }

    if (unknown)
        return GR_UNABLE;
    else
        return GR_DOMAIN;
}

int
gr_mat_find_nonzero_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    return GR_MAT_PIVOT_OP(ctx, MAT_FIND_NONZERO_PIVOT)(pivot_row, mat, start_row, end_row, column, ctx);
}
