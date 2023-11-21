/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: check that dimensions are compatible */
int
gr_mat_nonsingular_solve_lu_precomp(gr_mat_t X, const slong * perm,
    const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, c, n, m;

    n = gr_mat_nrows(X, ctx);
    m = gr_mat_ncols(X, ctx);

    if (X == B)
    {
        gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
        gr_ptr tmp = flint_malloc(sz * n);

        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
                set_shallow(GR_ENTRY(tmp, i, sz), GR_MAT_ENTRY(B, perm[i], c, sz), ctx);
            for (i = 0; i < n; i++)
                set_shallow(GR_MAT_ENTRY(X, i, c, sz), GR_ENTRY(tmp, i, sz), ctx);
        }

        flint_free(tmp);
    }
    else
    {
        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
            {
                status |= gr_set(GR_MAT_ENTRY(X, i, c, sz),
                    GR_MAT_ENTRY(B, perm[i], c, sz), ctx);
            }
        }
    }

    /* todo: inline for small n? */
    status |= gr_mat_nonsingular_solve_tril(X, A, X, 1, ctx);
    status |= gr_mat_nonsingular_solve_triu(X, A, X, 0, ctx);

    return status;
}
