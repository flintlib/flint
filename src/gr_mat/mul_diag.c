/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_mul_diag(gr_mat_t C, const gr_mat_t A, const gr_vec_t D, gr_ctx_t ctx)
{
    slong ar, ac, i, sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_srcptr Dptr = D->entries;

    ar = gr_mat_nrows(A, ctx);
    ac = gr_mat_ncols(A, ctx);

    if (ar != D->length || ar != gr_mat_nrows(C, ctx) || ac != gr_mat_ncols(C, ctx))
        return GR_DOMAIN;

    for (i = 0; i < ar; i++)
        status |= _gr_vec_mul(GR_MAT_ENTRY(C, i, 0, sz), GR_MAT_ENTRY(A, i, 0, sz), Dptr, ac, ctx);

    return status;

}
