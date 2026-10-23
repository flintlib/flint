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

truth_t
ca_mat_nonsingular_solve_lu(ca_mat_t X, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)
{
    int status;
    gr_ctx_t gr_ctx;
    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);

    status = gr_mat_nonsingular_solve_lu((gr_mat_struct *) X,
        (const gr_mat_struct *) A, (const gr_mat_struct *) B, gr_ctx);

    if (status & GR_UNABLE)
        return T_UNKNOWN;
    else if (status & GR_DOMAIN)
        return T_FALSE;
    else
        return T_TRUE;
}
