/*
    Copyright (C) 2011, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"

slong
nmod_mat_lu_with_pivots(slong * P, slong * pivots_nonpivots, nmod_mat_t A)
{
    slong i, j, k, n, rank;
    slong * pivots;
    slong * nonpivots;

    n = A->c;

    rank = nmod_mat_lu(P, A, 0);

    if (rank == 0)
    {
        for (i = 0; i < n; i++)
            pivots_nonpivots[i] = i;
        return rank;
    }

    pivots = pivots_nonpivots;
    nonpivots = pivots_nonpivots + rank;

    for (i = j = k = 0; i < rank; i++)
    {
        while (nmod_mat_entry(A, i, j) == UWORD(0))
        {
            nonpivots[k] = j;
            k++;
            j++;
        }
        pivots[i] = j;
        j++;
    }

    while (k < n - rank)
    {
        nonpivots[k] = j;
        k++;
        j++;
    }

    return rank;
}

