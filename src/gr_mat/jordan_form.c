/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_jordan_form(gr_mat_t J, gr_mat_t P, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_vec_t lambda;
    slong n;
    slong num_blocks, * block_size, * block_lambda;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(A, ctx);

    if (J == A || P == A)
    {
        gr_mat_t T;
        gr_mat_init(T, n, n, ctx);
        status |= gr_mat_set(T, A, ctx);
        status |= gr_mat_jordan_form(J, P, T, ctx);
        gr_mat_clear(T, ctx);
        return status;
    }

    gr_vec_init(lambda, 0, ctx);
    block_lambda = flint_malloc(sizeof(slong) * n);
    block_size = flint_malloc(sizeof(slong) * n);

    status = gr_mat_jordan_blocks(lambda, &num_blocks, block_lambda, block_size, A, ctx);

    if ((status == GR_SUCCESS) && P != NULL)
        status |= gr_mat_jordan_transformation(P, lambda, num_blocks, block_lambda, block_size, A, ctx);

    if (status == GR_SUCCESS)
        status |= gr_mat_set_jordan_blocks(J, lambda, num_blocks, block_lambda, block_size, ctx);

    gr_vec_clear(lambda, ctx);
    flint_free(block_lambda);
    flint_free(block_size);

    return status;
}
