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

#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_mat.h"

slong fmpz_mod_mat_nullspace(fmpz_mod_mat_t X, const fmpz_mod_mat_t A, const fmpz_mod_ctx_t ctx)
{
    slong i, j, k, m, n, rank, nullity;
    slong *p;
    slong *pivots;
    slong *nonpivots;
    fmpz_mod_mat_t tmp;

    m = A->r;
    n = A->c;

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    fmpz_mod_mat_init_set(tmp, A, ctx);
    rank = fmpz_mod_mat_rref(tmp, tmp, ctx);
    nullity = n - rank;

    fmpz_mod_mat_zero(X, ctx);

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
                fmpz_mod_neg(fmpz_mod_mat_entry(X, pivots[j], i),
                             fmpz_mod_mat_entry(tmp, j, nonpivots[i]), ctx);
            }

            fmpz_one(fmpz_mod_mat_entry(X, nonpivots[i], i));
        }
    }

    flint_free(p);
    fmpz_mod_mat_clear(tmp, ctx);

    return nullity;
}
