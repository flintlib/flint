/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int
ca_mat_exp(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)
{
    int success;
    slong i, j, len, n;
    ca_mat_t P, Q, J;
    ca_vec_t lambda, f_lambda;
    slong num_blocks, num_lambda, offset;
    slong * block_lambda, * block_size;

    n = ca_mat_nrows(A);

    if (n != ca_mat_ncols(A))
        return 0;

    if (n == 0)
        return 1;

    success = 1;

    ca_mat_init(P, n, n, ctx);
    ca_mat_init(Q, n, n, ctx);
    ca_mat_init(J, n, n, ctx);
    block_lambda = flint_malloc(sizeof(slong) * n);
    block_size = flint_malloc(sizeof(slong) * n);
    ca_vec_init(lambda, 0, ctx);
    ca_vec_init(f_lambda, 0, ctx);

    success = ca_mat_jordan_blocks(lambda, &num_blocks, block_lambda, block_size, A, ctx);
    if (!success)
        goto cleanup;

    num_lambda = ca_vec_length(lambda, ctx);

    success = ca_mat_jordan_transformation(P, lambda, num_blocks, block_lambda, block_size, A, ctx);
    if (!success)
        goto cleanup;

    if (ca_mat_inv(Q, P, ctx) != T_TRUE)
    {
        success = 0;
        goto cleanup;
    }

    /* Evaluate Jordan blocks and build matrix */
    ca_vec_set_length(f_lambda, num_lambda, ctx);
    for (i = 0; i < num_lambda; i++)
        ca_exp(ca_vec_entry(f_lambda, i), ca_vec_entry(lambda, i), ctx);

    offset = 0;
    for (i = 0; i < num_blocks; i++)
    {
        len = block_size[i];

        ca_set(ca_mat_entry(J, offset, offset), ca_vec_entry(f_lambda, block_lambda[i]), ctx);

        if (len > 1)
        {
            for (j = 1; j < len; j++)
                ca_div_ui(ca_mat_entry(J, offset, offset + j), ca_mat_entry(J, offset, offset + j - 1), FLINT_MAX(j, 1), ctx);

            for (j = 1; j < len; j++)
                _ca_vec_set(ca_mat_entry(J, offset + j, offset + j), ca_mat_entry(J, offset + j - 1, offset + j - 1), len - j, ctx);
        }

        offset += block_size[i];
    }

    ca_mat_mul(res, P, J, ctx);
    ca_mat_mul(res, res, Q, ctx);

cleanup:
    ca_mat_clear(P, ctx);
    ca_mat_clear(Q, ctx);
    ca_mat_clear(J, ctx);
    ca_vec_clear(lambda, ctx);
    ca_vec_clear(f_lambda, ctx);
    flint_free(block_lambda);
    flint_free(block_size);

    return success;
}
