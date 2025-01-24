/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_mat.h"

int
nmod_mat_solve(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i, rank, *perm;
    nmod_mat_t LU;
    int result;

    if (A->r == 0 || B->c == 0)
        return 1;

    nmod_mat_init_set(LU, A);
    perm = flint_malloc(sizeof(slong) * A->r);
    for (i = 0; i < A->r; i++)
        perm[i] = i;

    rank = nmod_mat_lu(perm, LU, 1);

    if (rank == A->r)
    {
        nmod_mat_t PB;
        nmod_mat_init(PB, B->r, B->c, B->mod.n);

        for (i = 0; i < A->r; i++)
            _nmod_vec_set(nmod_mat_entry_ptr(PB, i, 0), nmod_mat_entry_ptr(B, perm[i], 0), B->c);

        nmod_mat_solve_tril(X, LU, PB, 1);
        nmod_mat_solve_triu(X, LU, X, 0);

        nmod_mat_clear(PB);
        result = 1;
    }
    else
    {
        result = 0;
    }

    nmod_mat_clear(LU);
    flint_free(perm);

    return result;
}
