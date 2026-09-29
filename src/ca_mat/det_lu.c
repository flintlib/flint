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

int
ca_mat_det_lu(ca_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);

    if (gr_mat_det_lu(res, (const gr_mat_struct *) A, gr_ctx) != GR_SUCCESS)
    {
        ca_unknown(res, ctx);
        return 0;
    }

    return 1;
}
