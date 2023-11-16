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
gr_mat_randrank(gr_mat_t mat, flint_rand_t state,
                           slong rank, gr_ctx_t ctx)
{
    slong i;
    gr_ptr diag;
    int status = GR_SUCCESS;
    int parity;

    if (rank < 0 || rank > mat->r || rank > mat->c)
        return GR_DOMAIN;

    GR_TMP_INIT_VEC(diag, rank, ctx);

    for (i = 0; i < rank; i++)
        status |= gr_randtest_not_zero(GR_ENTRY(diag, i, ctx->sizeof_elem), state, ctx);

    status |= gr_mat_randpermdiag(&parity, mat, state, diag, rank, ctx);

    GR_TMP_CLEAR_VEC(diag, rank, ctx);

    return status;
}
