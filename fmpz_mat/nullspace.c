/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"

slong
fmpz_mat_nullspace(fmpz_mat_t res, const fmpz_mat_t mat)
{
    slong i, j, k, n, rank, nullity;
    slong * pivots;
    slong * nonpivots;
    fmpz_mat_t tmp;
    fmpz_t den;

    n = mat->c;

    fmpz_mat_init_set(tmp, mat);
    fmpz_init(den);

    rank = fmpz_mat_rref(tmp, den, mat);
    nullity = n - rank;

    fmpz_mat_zero(res);
    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            fmpz_one(res->rows[i] + i);
    }
    else if (nullity)
    {
        pivots = flint_malloc(rank * sizeof(slong));
        nonpivots = flint_malloc(nullity * sizeof(slong));

        for (i = j = k = 0; i < rank; i++)
        {
            while (fmpz_is_zero(tmp->rows[i] + j))
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

        fmpz_set(den, tmp->rows[0] + pivots[0]);

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
                fmpz_set(res->rows[pivots[j]] + i, tmp->rows[j] + nonpivots[i]);
            fmpz_neg(res->rows[nonpivots[i]] + i, den);
        }

        flint_free(pivots);
        flint_free(nonpivots);
    }

    fmpz_clear(den);
    fmpz_mat_clear(tmp);

    return nullity;
}
