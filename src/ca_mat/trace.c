/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr.h"
#include "gr_mat.h"

void
ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);

    if (!ca_mat_is_square(mat))
    {
        flint_throw(FLINT_ERROR, "ca_mat_trace: a square matrix is required!\n");
    }

    GR_MUST_SUCCEED(gr_mat_trace(trace, (const gr_mat_struct *) mat, gr_ctx));
}
