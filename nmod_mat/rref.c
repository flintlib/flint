/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "perm.h"

slong
_nmod_mat_rref(nmod_mat_t A, slong * pivots_nonpivots, slong * P)
{
    slong i, j, k, n, rank;
    slong * pivots;
    slong * nonpivots;

    nmod_mat_t U, V;

    n = A->c;

    rank = nmod_mat_lu(P, A, 0);

    if (rank == 0)
    {
        for (i = 0; i < n; i++)
            pivots_nonpivots[i] = i;
        return rank;
    }

    /* Clear L */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            nmod_mat_entry(A, i, j) = UWORD(0);

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    nmod_mat_init(U, rank, rank, A->mod.n);
    nmod_mat_init(V, rank, n - rank, A->mod.n);

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

    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            nmod_mat_entry(U, j, i) = nmod_mat_entry(A, j, pivots[i]);
    }

    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            nmod_mat_entry(V, j, i) = nmod_mat_entry(A, j, nonpivots[i]);
    }

    nmod_mat_solve_triu(V, U, V, 0);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
            nmod_mat_entry(A, j, pivots[i]) = (i == j);
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
    {
        for (j = 0; j < rank; j++)
            nmod_mat_entry(A, j, nonpivots[i]) = nmod_mat_entry(V, j, i);
    }

    nmod_mat_clear(U);
    nmod_mat_clear(V);

    return rank;
}

slong
nmod_mat_rref(nmod_mat_t A)
{
    slong rank, * pivots_nonpivots, * P;

    if (nmod_mat_is_empty(A))
        return 0;

    if (A->r == 1)
    {
        mp_limb_t c, cinv;
        slong i, j;
        slong r = 0;

        for (i = 0; i < A->c; i++)
        {
            c = nmod_mat_entry(A, 0, i);
            if (c != 0)
            {
                r = 1;
                if (c == 1)
                    break;

                cinv = nmod_inv(c, A->mod);
                nmod_mat_set_entry(A, 0, i, 1);
                for (j = i + 1;j < A->c; j++)
                {
                    nmod_mat_set_entry(A, 0, j, nmod_mul(nmod_mat_get_entry(A, 0, j), cinv, A->mod));
                }
                break;
            }
        }
        return r;
    }

    pivots_nonpivots = flint_malloc(sizeof(slong) * A->c);
    P = _perm_init(nmod_mat_nrows(A));

    rank = _nmod_mat_rref(A, pivots_nonpivots, P);

    flint_free(pivots_nonpivots);
    _perm_clear(P);

    return rank;
}
