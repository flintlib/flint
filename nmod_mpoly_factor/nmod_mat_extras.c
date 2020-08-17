/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int nmod_mat_is_reduced(const nmod_mat_t N)
{
    slong i, j, k = 0;
    slong r = nmod_mat_ncols(N);
    slong d = nmod_mat_nrows(N);
    
    for (i = 0; i < d; i++)
    for (j = 0; j < r; j++)
    {
        if (nmod_mat_entry(N, i, j) != 0)
        {
            if (nmod_mat_entry(N, i, j) == 1)
                k++;
            else
                return 0;
        }   
    }
    return k == r;
}

void nmod_mat_init_nullspace_tr(nmod_mat_t X, nmod_mat_t tmp)
{
    slong i, j, k, m, n, rank, nullity;
    slong * p;
    slong * pivots;
    slong * nonpivots;

    m = tmp->r;
    n = tmp->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    rank = nmod_mat_rref(tmp);

    nullity = n - rank;

    nmod_mat_init(X, nullity, n, tmp->mod.n);

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
                nmod_mat_entry(X, i, pivots[j]) = nmod_neg(c, tmp->mod);
            }

            nmod_mat_entry(X, i, nonpivots[i]) = UWORD(1);
        }
    }

    flint_free(p);
}

