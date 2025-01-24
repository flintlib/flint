/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
#include "templates.h"

int
TEMPLATE(T, mat_can_solve)(TEMPLATE(T, mat_t) X, const TEMPLATE(T, mat_t) A,
                           const TEMPLATE(T, mat_t) B, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j, k, col, *pivots, rank, *perm;
    TEMPLATE(T, mat_t) LU, LU2, PB;
    int result = 1;

    if (X->r != A->c || X->c != B->c)
        return 0;

    if (A->r == 0 || B->c == 0)
    {
        TEMPLATE(T, mat_zero)(X, ctx);

        return 1;
    }

    if (A->c == 0)
    {
       TEMPLATE(T, mat_zero)(X, ctx);

       return TEMPLATE(T, mat_is_zero)(B, ctx);
    }

    TEMPLATE(T, mat_init_set)(LU, A, ctx);
    perm = flint_malloc(sizeof(slong) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = TEMPLATE(T, mat_lu)(perm, LU, 0, ctx);

    TEMPLATE(T, mat_init)(PB, B->r, B->c, ctx);
    for (i = 0; i < B->r; i++)
        _TEMPLATE(T, vec_set)(TEMPLATE(T, mat_entry)(PB, i, 0), TEMPLATE(T, mat_entry)(B, perm[i], 0), B->c, ctx);

    TEMPLATE(T, mat_init)(LU2, rank, rank, ctx);

    pivots = flint_malloc(sizeof(slong)*rank);

    col = 0;
    for (i = 0; i < rank; i++)
    {
        while (TEMPLATE(T, is_zero)(TEMPLATE(T, mat_entry)(LU, i, col), ctx))
            col++;

        pivots[i] = col;

        for (j = 0; j < rank; j++)
            TEMPLATE(T, mat_entry_set)(LU2, j, i, TEMPLATE(T, mat_entry)(LU, j, col), ctx);

        col++;
    }

    X->r = rank;
    PB->r = rank;
    LU->r = rank;
    TEMPLATE(T, mat_solve_tril)(X, LU, PB, 1, ctx);
    LU->r = A->r;

    if (A->r > rank)
    {
        TEMPLATE(T, mat_t) P;

        LU->entries += rank * LU->stride;
        LU->r = A->r - rank;
        X->r = LU->c;

        TEMPLATE(T, mat_init)(P, LU->r, B->c, ctx);

        TEMPLATE(T, mat_mul)(P, LU, X, ctx);

        PB->r = LU->r;
        PB->entries += rank * PB->stride;

        result = TEMPLATE(T, mat_equal)(P, PB, ctx);

        PB->entries -= rank * PB->stride;
        TEMPLATE(T, mat_clear)(P, ctx);

        LU->entries -= rank * LU->stride;

        if (!result)
        {
            X->r = A->c;
            TEMPLATE(T, mat_zero)(X, ctx);
            goto cleanup;
        }
    }

    TEMPLATE(T, mat_solve_triu)(X, LU2, X, 0, ctx);

    X->r = A->c;

    k = rank - 1;
    for (i = A->c - 1; i >= 0; i--)
    {
        if (k < 0 || i != pivots[k])
        {
            for (j = 0; j < B->c; j++)
                TEMPLATE(T, zero)(TEMPLATE(T, mat_entry)(X, i, j), ctx);
        } else
        {
            for (j = 0; j < B->c; j++)
                TEMPLATE(T, mat_entry_set)(X, i, j, TEMPLATE(T, mat_entry)(X, k, j), ctx);

            k--;
        }
    }

cleanup:

    TEMPLATE(T, mat_clear)(LU2, ctx);

    PB->r = B->r;
    TEMPLATE(T, mat_clear)(PB, ctx);

    LU->r = A->r;
    TEMPLATE(T, mat_clear)(LU, ctx);
    flint_free(perm);

    flint_free(pivots);

    return result;
}

#endif
