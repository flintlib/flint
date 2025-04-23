/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_randsimilar(gr_mat_t mat, flint_rand_t state, slong opcount, gr_ctx_t ctx)
{
    slong c, i, j;
    slong m = mat->r;
    slong n = mat->c;
    slong sz = ctx->sizeof_elem;
    gr_mat_t W1, W2;
    gr_ptr x;
    int status = GR_SUCCESS;

    if (n != m)
        return GR_DOMAIN;

    if (n <= 1)
        return GR_SUCCESS;

    GR_TMP_INIT(x, ctx);

    /* Conjugate by random elementary operations. */
    for (c = 0; c < opcount; c++)
    {
        switch (n_randint(state, 3)) {
            case 0: {
                /* Swap row i and row j, swap column i and column j. */
                if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                    continue;
                status |= gr_mat_swap_rows(mat, NULL, i, j, ctx);
                status |= gr_mat_swap_cols(mat, NULL, i, j, ctx);
            } break;
            case 1: {
                /* For a random invertible x, scale row i by x, scale column i by the inverse of x. */
                i = n_randint(state, n);
                status |= gr_randtest_invertible(x, state, ctx);
                status |= _gr_vec_mul_scalar(GR_MAT_ENTRY(mat, i, 0, sz), GR_MAT_ENTRY(mat, i, 0, sz), mat->c, x, ctx);
                gr_mat_window_init(W1, mat, 0, i, n, i + 1, ctx);
                status |= gr_inv(x, x, ctx);
                status |= gr_mat_mul_scalar(W1, W1, x, ctx);
            } break;
            case 2: {
                /* For a random x, add x times row j to row i, subtract x times column i from column j. */
                if ((i = n_randint(state, n)) == (j = n_randint(state, n)))
                    continue;
                status |= gr_randtest(x, state, ctx);
                status |= _gr_vec_addmul_scalar(GR_MAT_ENTRY(mat, i, 0, sz), GR_MAT_ENTRY(mat, j, 0, sz), mat->c, x, ctx);
                gr_mat_window_init(W1, mat, 0, j, n, j + 1, ctx);
                gr_mat_window_init(W2, mat, 0, i, n, i + 1, ctx);
                status |= gr_mat_submul_scalar(W1, W2, x, ctx);
            } break;
            default: {
                /* Cannot happen. */
            } break;
        }
    }

    GR_TMP_CLEAR(x, ctx);

    return status;
}
