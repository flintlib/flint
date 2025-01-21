/*
    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_can_solve_inner(slong * rank, slong * prm, slong * piv,
                           nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, j, k, col, * pivots, * perm, rnk;
    nmod_mat_t LU, LU2, PB;
    int result = 1;

    if (rank == NULL)
        rank = & rnk;

    if (A->r != B->r || A->c != X->r || X->c != B->c)
    {
        *rank = 0;
        return 0;
    }

    if (A->r == 0 || B->c == 0)
    {
        nmod_mat_zero(X);

        *rank = 0;

        return 1;
    }

    if (A->c == 0)
    {
        nmod_mat_zero(X);

        *rank = 0;

        return nmod_mat_is_zero(B);
    }

    nmod_mat_init_set(LU, A);
    if (prm == NULL)
    {
        perm = flint_malloc(sizeof(slong)*A->r);
        for (i = 0; i < A->r; i++)
            perm[i] = i;
    } else
        perm = prm;

    *rank = nmod_mat_lu(perm, LU, 0);

    nmod_mat_init(PB, B->r, B->c, A->mod.n);
    for (i = 0; i < B->r; i++)
        _nmod_vec_set(nmod_mat_entry_ptr(PB, i, 0), nmod_mat_entry_ptr(B, perm[i], 0), B->c);

    nmod_mat_init(LU2, *rank, *rank, A->mod.n);

    if (piv == NULL)
        pivots = flint_malloc(sizeof(slong)*(*rank));
    else
        pivots = piv;

    col = 0;
    for (i = 0; i < *rank; i++)
    {
       while (nmod_mat_entry(LU, i, col) == 0)
           col++;

       pivots[i] = col;

       for (j = 0; j < *rank; j++)
           nmod_mat_set_entry(LU2, j, i, nmod_mat_entry(LU, j, col));

       col++;
    }

    X->r = *rank;
    PB->r = *rank;
    LU->r = *rank;
    nmod_mat_solve_tril(X, LU, PB, 1);
    LU->r = A->r;

    if (A->r > *rank)
    {
        nmod_mat_t P;

        LU->entries += *rank * LU->stride;
        LU->r = A->r - *rank;
        X->r = LU->c;

        nmod_mat_init(P, LU->r, B->c, A->mod.n);

        nmod_mat_mul(P, LU, X);

        PB->r = LU->r;
        PB->entries += *rank * PB->stride;

        result = nmod_mat_equal(P, PB);

        PB->entries -= *rank * PB->stride;
        nmod_mat_clear(P);

        LU->entries -= *rank * LU->stride;

        if (!result)
        {
            X->r = A->c;
            nmod_mat_zero(X);
            goto cleanup;
        }
    }

    nmod_mat_solve_triu(X, LU2, X, 0);

    X->r = A->c;

    k = (*rank) - 1;
    for (i = A->c - 1; i >= 0; i--)
    {
        if (k < 0 || i != pivots[k])
        {
            for (j = 0; j < B->c; j++)
                nmod_mat_set_entry(X, i, j, 0);
        } else
        {
            for (j = 0; j < B->c; j++)
                nmod_mat_set_entry(X, i, j, nmod_mat_entry(X, k, j));

            k--;
        }
    }

cleanup:

    nmod_mat_clear(LU2);

    PB->r = B->r;
    nmod_mat_clear(PB);

    LU->r = A->r;
    nmod_mat_clear(LU);
    if (prm == NULL)
        flint_free(perm);

    if (piv == NULL)
        flint_free(pivots);

    return result;
}

int
nmod_mat_can_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
{
    return nmod_mat_can_solve_inner(NULL, NULL, NULL, X, A, B);
}
