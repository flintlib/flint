/*
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

int fmpz_mod_mat_can_solve(fmpz_mod_mat_t X, const fmpz_mod_mat_t A,
                           const fmpz_mod_mat_t B)
{
    slong i, j, k, col, *pivots, rank, *perm;
    fmpz_mod_mat_t LU, LU2, PB;
    int result = 1;
    if (A->mat->r == 0 || B->mat->c == 0)
    {
        fmpz_mod_mat_zero(X);
                
        return 1;
    }

    if (A->mat->c == 0)
    {
       fmpz_mod_mat_zero(X);

       return fmpz_mod_mat_is_zero(B);
    }

    fmpz_mod_mat_init_set(LU, A);
    perm = flint_malloc(sizeof(slong) * A->mat->r);
    for (i = 0; i < A->mat->r; i++)
        perm[i] = i;

    rank = fmpz_mod_mat_lu(perm, LU, 0);

    fmpz_mod_mat_window_init(PB, B, 0, 0, B->mat->r, B->mat->c);
    for (i = 0; i < B->mat->r; i++)
        PB->mat->rows[i] = B->mat->rows[perm[i]];

    fmpz_mod_mat_init(LU2, rank, rank, A->mod);

    pivots = flint_malloc(sizeof(slong)*rank);

    col = 0;
    for (i = 0; i < rank; i++)
    {
        while (fmpz_is_zero(fmpz_mod_mat_entry(LU, i, col)))
            col++;

        pivots[i] = col;

        for (j = 0; j < rank; j++)
            fmpz_mod_mat_set_entry(LU2, j, i, fmpz_mod_mat_entry(LU, j, col));

        col++;
    }

    X->mat->r = rank;
    PB->mat->r = rank;
    LU->mat->r = rank;
    fmpz_mod_mat_solve_tril(X, LU, PB, 1);
    LU->mat->r = A->mat->r;

    if (A->mat->r > rank)
    {
        fmpz_mod_mat_t P;

        LU->mat->rows += rank;
        LU->mat->r = A->mat->r - rank;

        fmpz_mod_mat_init(P, LU->mat->r, B->mat->c, A->mod);

        fmpz_mod_mat_mul(P, LU, X);

        PB->mat->r = LU->mat->r;
        PB->mat->rows += rank;

        result = fmpz_mod_mat_equal(P, PB);

        PB->mat->rows -= rank;
        fmpz_mod_mat_clear(P);

        LU->mat->rows -= rank;

        if (!result)
        {
            fmpz_mod_mat_zero(X);
            goto cleanup;
        }
    }

    fmpz_mod_mat_solve_triu(X, LU2, X, 0);

    X->mat->r = A->mat->c;

    k = rank - 1;
    for (i = A->mat->c - 1; i >= 0; i--)
    {
        if (k < 0 || i != pivots[k])
        {
            for (j = 0; j < B->mat->c; j++)
                fmpz_zero(fmpz_mod_mat_entry(X, i, j));
        }
        else
        {
            for (j = 0; j < B->mat->c; j++)
                fmpz_mod_mat_set_entry(X, i, j, fmpz_mod_mat_entry(X, k, j));

            k--;
        }
    }

cleanup:

    fmpz_mod_mat_clear(LU2);
    
    PB->mat->r = B->mat->r;
    fmpz_mod_mat_window_clear(PB);

    fmpz_mod_mat_clear(LU);
    flint_free(perm);

    flint_free(pivots);

    return result;
}

