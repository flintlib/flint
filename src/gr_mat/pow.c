/*
    Copyright (C) 2025 Lars Göttgens

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "gr.h"
#include "gr_mat.h"

int
gr_mat_pow_ui(gr_mat_t res, const gr_mat_t mat, ulong exp, gr_ctx_t ctx)
{
    int status;
    slong sz = ctx->sizeof_elem;
    slong d;

    d = gr_mat_nrows(res, ctx);

    if (d != gr_mat_ncols(res, ctx) || d != gr_mat_nrows(mat, ctx)
        || d != gr_mat_ncols(mat, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;

    if (exp <= 2 || d <= 1)
    {
        if (exp == 0 || d == 0)
        {
            status |= gr_mat_one(res, ctx);
        }
        else if (d == 1)
        {
            status |= gr_pow_ui(GR_MAT_ENTRY(res, 0, 0, sz),
                                GR_MAT_ENTRY(mat, 0, 0, sz), exp, ctx);
        }
        else if (exp == 1)
        {
            status |= gr_mat_set(res, mat, ctx);
        }
        else if (exp == 2)
        {
            status |= gr_mat_sqr(res, mat, ctx);
        }
    }
    else
    {
        gr_ptr T, U;
        slong i;

        status |= gr_mat_init_set(T, mat, ctx);
        gr_mat_init(U, d, d, ctx);

        for (i = FLINT_BIT_COUNT(exp) - 2; i >= 0; i--)
        {
            status |= gr_mat_sqr(U, T, ctx);

            if (exp & (WORD(1) << i))
                status |= gr_mat_mul(T, U, mat, ctx);
            else
                gr_mat_swap(T, U, ctx);
        }

        gr_mat_swap(res, T, ctx);
        gr_mat_clear(T, ctx);
        gr_mat_clear(U, ctx);
    }

    return status;
}
