/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_exp_jet(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
{
    slong i;
    int status;
    slong sz = ctx->sizeof_elem;

    if (len <= 0)
        return GR_SUCCESS;

    status = gr_exp(res, x, ctx);

    if (status == GR_SUCCESS)
    {
        for (i = 1; i < len; i++)
            status |= gr_div_ui(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), i, ctx);
    }

    return status;
}

int
gr_mat_exp_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    return gr_mat_func_jordan(res, A, (gr_method_vec_op) gr_exp_jet, ctx);
}

int
gr_mat_exp(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    return GR_MAT_UNARY_OP(ctx, MAT_EXP)(res, A, ctx);
}
