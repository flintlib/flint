/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"
#include "gr_mat.h"

truth_t
ca_mat_log(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    gr_ctx_t cctx;
    int status;
    _gr_ctx_init_ca_from_ref(cctx, GR_CTX_CC_CA, ctx);
    status = gr_mat_log((gr_mat_struct *) res, (const gr_mat_struct *) A, cctx);
    if (status & GR_UNABLE)
        return T_UNKNOWN;
    if (status & GR_DOMAIN)
        return T_FALSE;
    return T_TRUE;
}
