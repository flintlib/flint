/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart
    Copyright (C) 2022 Fredrik Johansson

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
gr_mat_solve_field(gr_mat_t X, const gr_mat_t A,
                           const gr_mat_t B, gr_ctx_t ctx)
{
    slong i, j, k, col, *pivots, rank, *perm;
    gr_mat_t LU, LU2, PB;
    int status = GR_SUCCESS;
    truth_t is_zero;
    slong sz = ctx->sizeof_elem;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        return GR_DOMAIN;
    }

    if (A->r == 0 || B->c == 0)
    {
        return gr_mat_zero(X, ctx);
    }

    if (A->c == 0)
    {
        status = gr_mat_zero(X, ctx);
        if (status != GR_SUCCESS)
            return status;

        is_zero = gr_mat_is_zero(B, ctx);
        if (is_zero == T_TRUE)
            return GR_SUCCESS;
        if (is_zero == T_FALSE)
            return GR_DOMAIN;
        return GR_UNABLE;
    }

    status |= gr_mat_init_set(LU, A, ctx);

    perm = flint_malloc(sizeof(slong) * A->r);

    for (i = 0; i < A->r; i++)
        perm[i] = i;

    status |= gr_mat_lu(&rank, perm, LU, LU, 0, ctx);
    if (status != GR_SUCCESS)
        goto cleanup1;

    gr_mat_init(PB, B->r, B->c, ctx);
    for (i = 0; i < B->r; i++)
        status |= _gr_vec_set(GR_MAT_ENTRY(PB, i, 0, sz), GR_MAT_ENTRY(B, perm[i], 0, sz), B->c, ctx);

    gr_mat_init(LU2, rank, rank, ctx);
    pivots = flint_malloc(sizeof(slong) * rank);

    col = 0;
    for (i = 0; i < rank; i++)
    {
        while (1)
        {
            is_zero = gr_is_zero(gr_mat_entry_ptr(LU, i, col, ctx), ctx);

            if (is_zero == T_UNKNOWN)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            if (is_zero == T_TRUE)
            {
                col++;
            }
            else
            {
                break;
            }
        }

        pivots[i] = col;

        for (j = 0; j < rank; j++)
            status |= gr_set(gr_mat_entry_ptr(LU2, j, i, ctx), gr_mat_entry_ptr(LU, j, col, ctx), ctx);

        col++;
    }

    X->r = rank;
    PB->r = rank;
    LU->r = rank;
    status |= gr_mat_nonsingular_solve_tril(X, LU, PB, 1, ctx);
    if (status != GR_SUCCESS)
        goto cleanup;

    LU->r = A->r;

    if (A->r > rank)
    {
        gr_mat_t P;
        truth_t equal;

        /* LU->rows += rank */
        LU->entries = GR_ENTRY(LU->entries, rank * LU->stride, sz);
        LU->r = A->r - rank;
        X->r = LU->c;

        gr_mat_init(P, LU->r, B->c, ctx);

        status |= gr_mat_mul(P, LU, X, ctx);

        PB->r = LU->r;
        /* PB->rows += rank */
        PB->entries = GR_ENTRY(PB->entries, rank * PB->stride, sz);

        equal = gr_mat_equal(P, PB, ctx);

        /* PB->rows -= rank; */
        PB->entries = GR_ENTRY(PB->entries, -rank * PB->stride, sz);
        gr_mat_clear(P, ctx);

        /* LU->rows -= rank; */
        LU->entries = GR_ENTRY(LU->entries, -rank * LU->stride, sz);
        LU->r = A->r;

        if (equal == T_UNKNOWN)
        {
            X->r = A->c;
            status |= GR_UNABLE;
            goto cleanup;
        }

        if (status == GR_SUCCESS && equal == T_FALSE)
        {
            /* restore original size of X */
            X->r = A->c;
            status = GR_DOMAIN;
            status |= gr_mat_zero(X, ctx);
            goto cleanup;
        }
    }

    status |= gr_mat_nonsingular_solve_triu(X, LU2, X, 0, ctx);
    if (status != GR_SUCCESS)
    {
        X->r = A->c;
        goto cleanup;
    }

    /* restore original size of X */
    X->r = A->c;

    k = rank - 1;
    for (i = A->c - 1; i >= 0; i--)
    {
        if (k < 0 || i != pivots[k])
        {
            for (j = 0; j < B->c; j++)
                status |= gr_zero(gr_mat_entry_ptr(X, i, j, ctx), ctx);
        }
        else
        {
            for (j = 0; j < B->c; j++)
                status |= gr_set(gr_mat_entry_ptr(X, i, j, ctx), gr_mat_entry_ptr(X, k, j, ctx), ctx);

            k--;
        }
    }

cleanup:

    gr_mat_clear(LU2, ctx);

    PB->r = B->r;
    gr_mat_clear(PB, ctx);

    flint_free(pivots);

cleanup1:
    /* restore original size of LU */
    LU->r = A->r;
    gr_mat_clear(LU, ctx);
    flint_free(perm);

    return status;
}
