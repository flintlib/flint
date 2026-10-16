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
ca_mat_transpose(ca_mat_t B, const ca_mat_t A, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);

    if (ca_mat_nrows(B) != ca_mat_ncols(A) || ca_mat_ncols(B) != ca_mat_nrows(A))
    {
        flint_throw(FLINT_ERROR, "Exception (ca_mat_transpose). Incompatible dimensions.\n");
    }

    GR_MUST_SUCCEED(gr_mat_transpose((gr_mat_struct *) B, (const gr_mat_struct *) A, gr_ctx));
}
