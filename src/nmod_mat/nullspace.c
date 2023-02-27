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

slong
nmod_mat_nullspace(nmod_mat_t X, const nmod_mat_t A)
{
    slong i, j, k, m, n, rank, nullity;
    slong * p;
    slong * pivots;
    slong * nonpivots;
    nmod_mat_t tmp;

    m = A->r;
    n = A->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    nmod_mat_init_set(tmp, A);
    rank = nmod_mat_rref(tmp);
    nullity = n - rank;

    nmod_mat_zero(X);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            nmod_mat_entry(X, i, i) = UWORD(1);
    }
    else if (nullity)
    {
        pivots = p;            /* length = rank */
        nonpivots = p + rank;  /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (nmod_mat_entry(tmp, i, j) == UWORD(0))
            {
                nonpivots[k] = j;
                k++;
                j++;
            }
            pivots[i] = j;
            j++;
        }
        while (k < nullity)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
            {
                mp_limb_t c = nmod_mat_entry(tmp, j, nonpivots[i]);
                nmod_mat_entry(X, pivots[j], i) = nmod_neg(c, A->mod);
            }

            nmod_mat_entry(X, nonpivots[i], i) = UWORD(1);
        }
    }

    flint_free(p);
    nmod_mat_clear(tmp);

    return nullity;
}
