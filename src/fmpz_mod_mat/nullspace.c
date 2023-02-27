/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mat.h"

slong fmpz_mod_mat_nullspace(fmpz_mod_mat_t X, const fmpz_mod_mat_t A)
{
    slong i, j, k, m, n, rank, nullity;
    slong *p;
    slong *pivots;
    slong *nonpivots;
    fmpz_mod_mat_t tmp;

    m = A->mat->r;
    n = A->mat->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    fmpz_mod_mat_init_set(tmp, A);
    rank = fmpz_mod_mat_rref(NULL, tmp);
    nullity = n - rank;

    fmpz_mod_mat_zero(X);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            fmpz_one(fmpz_mod_mat_entry(X, i, i));
    }
    else if (nullity)
    {
        pivots = p;             /* length = rank */
        nonpivots = p + rank;   /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (fmpz_is_zero(fmpz_mod_mat_entry(tmp, i, j)))
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
                fmpz_negmod(fmpz_mod_mat_entry(X, pivots[j], i),
                             fmpz_mod_mat_entry(tmp, j, nonpivots[i]), A->mod);
            }

            fmpz_one(fmpz_mod_mat_entry(X, nonpivots[i], i));
        }
    }

    flint_free(p);
    fmpz_mod_mat_clear(tmp);

    return nullity;
}

