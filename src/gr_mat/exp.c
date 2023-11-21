/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_special.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_exp_jordan(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, len, n;
    gr_mat_t P, Q, J;
    gr_vec_t lambda, f_lambda;
    slong num_blocks, num_lambda, offset;
    slong * block_lambda, * block_size;

    n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return GR_DOMAIN;

    if (n == 0)
        return GR_SUCCESS;

    gr_mat_init(P, n, n, ctx);
    gr_mat_init(Q, n, n, ctx);
    gr_mat_init(J, n, n, ctx);

    block_lambda = flint_malloc(sizeof(slong) * n);
    block_size = flint_malloc(sizeof(slong) * n);
    gr_vec_init(lambda, 0, ctx);
    gr_vec_init(f_lambda, 0, ctx);

    status |= gr_mat_jordan_blocks(lambda, &num_blocks, block_lambda, block_size, A, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    num_lambda = gr_vec_length(lambda, ctx);

    status |= gr_mat_jordan_transformation(P, lambda, num_blocks, block_lambda, block_size, A, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    status |= gr_mat_inv(Q, P, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    /* Evaluate Jordan blocks and build matrix */
    gr_vec_set_length(f_lambda, num_lambda, ctx);
    for (i = 0; i < num_lambda && status == GR_SUCCESS; i++)
        status |= gr_exp(gr_vec_entry_ptr(f_lambda, i, ctx), gr_vec_entry_srcptr(lambda, i, ctx), ctx);

    offset = 0;
    for (i = 0; i < num_blocks; i++)
    {
        len = block_size[i];

        status |= gr_set(gr_mat_entry_ptr(J, offset, offset, ctx), gr_vec_entry_srcptr(f_lambda, block_lambda[i], ctx), ctx);

        if (len > 1)
        {
            for (j = 1; j < len; j++)
                status |= gr_div_ui(gr_mat_entry_ptr(J, offset, offset + j, ctx), gr_mat_entry_srcptr(J, offset, offset + j - 1, ctx), FLINT_MAX(j, 1), ctx);

            for (j = 1; j < len; j++)
                status |= _gr_vec_set(gr_mat_entry_ptr(J, offset + j, offset + j, ctx), gr_mat_entry_srcptr(J, offset + j - 1, offset + j - 1, ctx), len - j, ctx);
        }

        offset += block_size[i];
    }

    status |= gr_mat_mul(res, P, J, ctx);
    status |= gr_mat_mul(res, res, Q, ctx);

cleanup:
    gr_mat_clear(P, ctx);
    gr_mat_clear(Q, ctx);
    gr_mat_clear(J, ctx);
    gr_vec_clear(lambda, ctx);
    gr_vec_clear(f_lambda, ctx);
    flint_free(block_lambda);
    flint_free(block_size);

    return status;
}

int
gr_mat_exp(gr_mat_t res, const gr_mat_t A, gr_ctx_t ctx)
{
    return GR_MAT_UNARY_OP(ctx, MAT_EXP)(res, A, ctx);
}
