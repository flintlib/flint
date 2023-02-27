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
ca_mat_jordan_form(ca_mat_t J, ca_mat_t P, const ca_mat_t A, ca_ctx_t ctx)
{
    ca_vec_t lambda;
    slong n;
    slong num_blocks, * block_size, * block_lambda;
    int success;

    n = ca_mat_nrows(A);

    if (J == A || P == A)
    {
        ca_mat_t T;
        ca_mat_init(T, n, n, ctx);
        ca_mat_set(T, A, ctx);
        success = ca_mat_jordan_form(J, P, T, ctx);
        ca_mat_clear(T, ctx);
        return success;
    }

    ca_vec_init(lambda, 0, ctx);
    block_lambda = flint_malloc(sizeof(slong) * n);
    block_size = flint_malloc(sizeof(slong) * n);

    success = ca_mat_jordan_blocks(lambda, &num_blocks, block_lambda, block_size, A, ctx);

    if (success && P != NULL)
        success = ca_mat_jordan_transformation(P, lambda, num_blocks, block_lambda, block_size, A, ctx);

    if (success)
        ca_mat_set_jordan_blocks(J, lambda, num_blocks, block_lambda, block_size, ctx);

    ca_vec_clear(lambda, ctx);
    flint_free(block_lambda);
    flint_free(block_size);

    return success;
}
