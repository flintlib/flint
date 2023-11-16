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

int
ca_mat_jordan_transformation(ca_mat_t mat, const ca_vec_t lambda, slong num_blocks, slong * block_lambda, slong * block_size, const ca_mat_t A, ca_ctx_t ctx)
{
    gr_ctx_t cctx;
    int success;
    _gr_ctx_init_ca_from_ref(cctx, GR_CTX_CC_CA, ctx);
    success = (gr_mat_jordan_transformation((gr_mat_struct *) mat, (gr_vec_struct *) lambda, num_blocks, block_lambda, block_size, (const gr_mat_struct *) A, cctx) == GR_SUCCESS);
    return success;
}
