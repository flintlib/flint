/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod_mat.h"

int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                           const fmpz_mod_mat_t B, const fmpz_mod_ctx_t ctx)
{
    slong i, j, k, col, *pivots, rank, *perm;
    fmpz_mod_mat_t LU, LU2, PB;
    int result = 1;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        return 0;
    }

    if (A->r == 0 || B->c == 0)
    {
        fmpz_mod_mat_zero(X, ctx);
        return 1;
    }

    if (A->c == 0)
    {
       fmpz_mod_mat_zero(X, ctx);
       return fmpz_mod_mat_is_zero(B, ctx);
    }

    fmpz_mod_mat_init_set(LU, A, ctx);
    perm = flint_malloc(sizeof(slong) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = fmpz_mod_mat_lu(perm, LU, 0, ctx);

    fmpz_mod_mat_window_init(PB, B, 0, 0, B->r, B->c, ctx);
    for (i = 0; i < B->r; i++)
        PB->rows[i] = B->rows[perm[i]];

    fmpz_mod_mat_init(LU2, rank, rank, ctx);

    pivots = flint_malloc(sizeof(slong)*rank);

    col = 0;
    for (i = 0; i < rank; i++)
    {
        while (fmpz_is_zero(fmpz_mod_mat_entry(LU, i, col)))
            col++;

        pivots[i] = col;

        for (j = 0; j < rank; j++)
            fmpz_mod_mat_set_entry(LU2, j, i, fmpz_mod_mat_entry(LU, j, col), ctx);

        col++;
    }

    X->r = rank;
    PB->r = rank;
    LU->r = rank;
    fmpz_mod_mat_solve_tril(X, LU, PB, 1, ctx);
    LU->r = A->r;

    if (A->r > rank)
    {
        fmpz_mod_mat_t P;

        LU->rows += rank;
        LU->r = A->r - rank;
        X->r = LU->c;

        fmpz_mod_mat_init(P, LU->r, B->c, ctx);

        fmpz_mod_mat_mul(P, LU, X, ctx);

        PB->r = LU->r;
        PB->rows += rank;

        result = fmpz_mod_mat_equal(P, PB, ctx);

        PB->rows -= rank;
        fmpz_mod_mat_clear(P, ctx);

        LU->rows -= rank;

        if (!result)
        {
            X->r = A->c;
            fmpz_mod_mat_zero(X, ctx);
            goto cleanup;
        }
    }

    fmpz_mod_mat_solve_triu(X, LU2, X, 0, ctx);

    X->r = A->c;

    k = rank - 1;
    for (i = A->c - 1; i >= 0; i--)
    {
        if (k < 0 || i != pivots[k])
        {
            for (j = 0; j < B->c; j++)
                fmpz_zero(fmpz_mod_mat_entry(X, i, j));
        }
        else
        {
            for (j = 0; j < B->c; j++)
                fmpz_mod_mat_set_entry(X, i, j, fmpz_mod_mat_entry(X, k, j), ctx);

            k--;
        }
    }

cleanup:

    fmpz_mod_mat_clear(LU2, ctx);

    PB->r = B->r;
    fmpz_mod_mat_window_clear(PB, ctx);

    LU->r = A->r;
    fmpz_mod_mat_clear(LU, ctx);
    flint_free(perm);

    flint_free(pivots);

    return result;
}
