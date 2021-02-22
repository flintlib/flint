/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


int fmpz_mod_mat_is_reduced(const fmpz_mod_mat_t N)
{
    slong i, j, k = 0;
    slong r = fmpz_mod_mat_ncols(N);
    slong d = fmpz_mod_mat_nrows(N);
    
    for (i = 0; i < d; i++)
    for (j = 0; j < r; j++)
    {
        if (!fmpz_is_zero(fmpz_mod_mat_entry(N, i, j)))
        {
            if (fmpz_is_one(fmpz_mod_mat_entry(N, i, j)))
                k++;
            else
                return 0;
        }   
    }
    return k == r;
}

void fmpz_mod_mat_init_nullspace_tr(fmpz_mod_mat_t X, fmpz_mod_mat_t tmp, const fmpz_mod_ctx_t ctx)
{
    slong i, j, k, m, n, rank, nullity;
    slong * p;
    slong * pivots;
    slong * nonpivots;

    m = fmpz_mod_mat_nrows(tmp);
    n = fmpz_mod_mat_ncols(tmp);

    p = FLINT_ARRAY_ALLOC(FLINT_MAX(m, n), slong);

    rank = fmpz_mod_mat_rref(NULL, tmp);

    nullity = n - rank;

    fmpz_mod_mat_init(X, nullity, n, fmpz_mod_ctx_modulus(ctx));

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            fmpz_one(fmpz_mod_mat_entry(X, i, i));
    }
    else if (nullity)
    {
        pivots = p;            /* length = rank */
        nonpivots = p + rank;  /* length = nullity */

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
                fmpz_mod_neg(fmpz_mod_mat_entry(X, i, pivots[j]),
                             fmpz_mod_mat_entry(tmp, j, nonpivots[i]), ctx);
            }

            fmpz_one(fmpz_mod_mat_entry(X, i, nonpivots[i]));
        }
    }

    flint_free(p);
}

