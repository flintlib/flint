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
ca_mat_nonsingular_lu(slong * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx)
{
    gr_ctx_t gr_ctx;
    slong rank;
    int status;

    if (ca_mat_is_empty(A))
        return T_TRUE;

    _gr_ctx_init_ca_from_ref(gr_ctx, GR_CTX_CC_CA, ctx);
    status = gr_mat_lu_generic(&rank, P, (gr_mat_struct *) LU,
        (const gr_mat_struct *) A, 1, gr_ctx);

    if (status != GR_SUCCESS)
        return T_UNKNOWN;

    if (rank == 0)
        return T_FALSE;

    return T_TRUE;
}
