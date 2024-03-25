/*
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2024 Kartik Venkatram

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "stdlib.h"
#include "gr_sparse_mat.h"

int gr_lil_mat_nullspace(gr_mat_t X, const gr_lil_mat_t M, flint_rand_t state, slong max_iters, const char *algorithm, slong block_size, gr_ctx_t ctx)
{
    /* Generate random solutions to a random system Mx = b and stop when nullspace filled */
    slong i, j, c, iter, nullity, *x_pivots, sz;
    gr_ptr x, coeff, *xs;
    int status = GR_SUCCESS, cur_status;
    
    sz = ctx->sizeof_elem;
    c = gr_sparse_mat_ncols(M, ctx);

    GR_TMP_INIT(coeff, ctx);
    GR_TMP_INIT_VEC(x, c, ctx);
    nullity = 0;
    xs = NULL;
    x_pivots = NULL;
    for (iter = 0; iter < max_iters; )
    {
        if (strcmp(algorithm, "lanczos") == 0)
            cur_status = gr_lil_mat_nullvector_lanczos(x, M, state, ctx);
        else if (strcmp(algorithm, "wiedemann") == 0)
            cur_status = gr_lil_mat_nullvector_wiedemann(x, M, state, ctx);
        else if (strcmp(algorithm, "block lanczos") == 0)
            cur_status = gr_lil_mat_nullvector_block_lanczos(x, M, block_size, state, ctx);
        else if (strcmp(algorithm, "block wiedemann") == 0)
            cur_status = gr_lil_mat_nullvector_block_wiedemann(x, M, block_size, state, ctx);
        else
        {
            status = GR_DOMAIN;
            break;
        }
        if (cur_status == GR_TEST_FAIL)
        {
            ++iter; 
            continue;
        }
        
        /* Reduce by existing kernel vectors */
        for (j = nullity-1; j >= 0; --j) 
        {
            status |= gr_neg(coeff, GR_ENTRY(x, x_pivots[j], sz), ctx);
            status |= _gr_vec_addmul_scalar(x, xs[j], c, coeff, ctx);
        }

        /* Normalize last nonzero entry to 1 */
        for (i = c-1; i >= 0; --i)
            if (gr_is_zero(GR_ENTRY(x, i, sz), ctx) != T_TRUE)
                break;
        if (i == -1) {
            /* x in span of xs, nullspace probably complete */
            ++iter;
            continue;
        }
        status |= gr_inv(coeff, GR_ENTRY(x, i, sz), ctx);
        status |= _gr_vec_mul_scalar(x, x, c, coeff, ctx);

        /* Reduce previous vectors by this one */
        for (j = 0; j < nullity; ++j) 
        {
            status |= gr_neg(coeff, GR_ENTRY(xs[j], i, sz), ctx);
            status |= _gr_vec_addmul_scalar(xs[j], x, c, coeff, ctx);
        }

        /* Insert into list of vectors in nullspace (ordered by pivot) */
        xs = flint_realloc(xs, (nullity+1)*sizeof(gr_ptr));
        x_pivots = flint_realloc(x_pivots, (nullity+1)*sizeof(slong *));
        for (j = 0; j < nullity; ++j)
            if (i > x_pivots[j])
                break;
        memmove(xs + j + 1, xs + j, (nullity - j) * sizeof(gr_ptr));
        memmove(x_pivots + j + 1, x_pivots + j, (nullity - j) * sizeof(slong *));
        x_pivots[j] = i;
        xs[j] = x; // Steal allocation
        x = NULL;
        GR_TMP_INIT_VEC(x, c, ctx);

        // Advance nullity and restart iteration
        nullity += 1;
        iter = 0;
    }
    // Set X to have xs as column vectors
    // TODO: can use shallow to reuse memory?
    gr_mat_init(X, c, nullity, ctx);
    for (i = 0; i < nullity; ++i)
        for (j = 0; j < c; ++j) 
            status |= gr_set(gr_mat_entry_ptr(X, j, i, ctx), GR_ENTRY(xs[i], j, sz), ctx);

    flint_free(x_pivots);
    GR_TMP_CLEAR_VEC(x, c, ctx);
    for (i = 0; i < nullity; ++i)
        GR_TMP_CLEAR_VEC(xs[i], c, ctx);
    flint_free(xs);
    return status;
}

