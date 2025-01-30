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
_gr_mat_func_jordan(gr_mat_t res, const gr_mat_t A, gr_method_vec_op jet_func1, gr_method_vec_scalar_op jet_func2, gr_srcptr func2_param, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, k, len, n;
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

    /* Over inexact implementations of R and C, the Jordan transformation
       has virtually no chance of succeeding, so we fall back on
       looking for a diagonalization with distinct eigenvalues.

        TODO: do the same in gr_mat_jordan_form
        TODO: can we do better than this?
    */
    if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        status = gr_mat_diagonalization(lambda, Q, P, A, 0, ctx);

        if (status != GR_SUCCESS)
            goto cleanup;

        num_blocks = num_lambda = n;
        for (i = 0; i < n; i++)
        {
            block_lambda[i] = i;
            block_size[i] = 1;
        }
    }
    else
    {
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
    }

    /* Evaluate Jordan blocks and build matrix */

    /* Sanity check to make sure that we can loop over the lambdas
       and then over the blocks in an inner loop and get the right
       offsets. */
    for (i = 1; i < num_blocks; i++)
    {
        if (block_lambda[i] < block_lambda[i - 1])
        {
            flint_printf("ERROR: jordan form: block indices not sorted\n");
            flint_abort();
        }
    }

    offset = 0;
    for (i = 0; i < num_lambda && status == GR_SUCCESS; i++)
    {
        len = 0;
        for (k = 0; k < num_blocks; k++)
        {
            if (i == block_lambda[k])
                len = FLINT_MAX(len, block_size[k]);
        }

        gr_vec_fit_length(f_lambda, len, ctx);

        if (jet_func1 != NULL)
            status |= jet_func1(f_lambda->entries, gr_vec_entry_srcptr(lambda, i, ctx), len, ctx);
        else
            status |= jet_func2(f_lambda->entries, gr_vec_entry_srcptr(lambda, i, ctx), len, func2_param, ctx);

        for (k = 0; k < num_blocks; k++)
        {
            if (i == block_lambda[k])
            {
                /* Copy shifted rows */
                len = block_size[k];
                for (j = 0; j < len; j++)
                    status |= _gr_vec_set(gr_mat_entry_ptr(J, offset + j, offset + j, ctx), f_lambda->entries, len - j, ctx);
                offset += block_size[k];
            }
        }
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
gr_mat_func_jordan(gr_mat_t res, const gr_mat_t A, gr_method_vec_op jet_func, gr_ctx_t ctx)
{
    return _gr_mat_func_jordan(res, A, jet_func, NULL, NULL, ctx);
}

int
gr_mat_func_param_jordan(gr_mat_t res, const gr_mat_t A, gr_method_vec_scalar_op jet_func, gr_srcptr c, gr_ctx_t ctx)
{
    return _gr_mat_func_jordan(res, A, NULL, jet_func, c, ctx);
}
