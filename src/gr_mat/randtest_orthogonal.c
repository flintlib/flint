/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
gr_mat_randtest_orthogonal(gr_mat_t A, flint_rand_t state, gr_ctx_t ctx)
{
    slong n = A->r;
    slong i, j, density;
    int status = GR_SUCCESS;

    if (n != A->c)
        return GR_DOMAIN;

    /* See https://arxiv.org/abs/math/0606320 */
    /* If S is skew-symmetric, then (I + S)^(-1) (I - S) is orthogonal. */
    gr_mat_t S, T;

    gr_mat_init(S, n, n, ctx);
    gr_mat_init(T, n, n, ctx);

    /* Generate a random skew-symmetric matrix */
    density = n_randint(state, 16);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (n_randint(state, 16) < density)
            {
                status |= gr_randtest(gr_mat_entry_ptr(S, i, j, ctx), state, ctx);
                status |= gr_neg(gr_mat_entry_ptr(S, j, i, ctx), gr_mat_entry_ptr(S, i, j, ctx), ctx);
            }
        }
    }

    status |= gr_mat_neg(T, S, ctx);
    status |= gr_mat_add_ui(S, S, 1, ctx);
    status |= gr_mat_add_ui(T, T, 1, ctx);
    status |= gr_mat_nonsingular_solve(A, S, T, ctx);

    gr_mat_clear(S, ctx);
    gr_mat_clear(T, ctx);

    /* Fall back to a random permutation matrix */
    if (status != GR_SUCCESS)
    {
        slong * P;
        status = gr_mat_zero(A, ctx);

        P = _perm_init(n);
        _perm_randtest(P, n, state);
        for (i = 0; i < n; i++)
            status |= gr_one(gr_mat_entry_ptr(A, i, P[i], ctx), ctx);
        _perm_clear(P);
    }

    for (i = 0; i < n; i++)
        if (n_randint(state, 2))
            status |= _gr_vec_neg(gr_mat_entry_ptr(A, i, 0, ctx),
                                  gr_mat_entry_ptr(A, i, 0, ctx), n, ctx);

    return status;
}

