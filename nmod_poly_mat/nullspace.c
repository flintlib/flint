/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"

slong
nmod_poly_mat_nullspace(nmod_poly_mat_t res, const nmod_poly_mat_t mat)
{
    slong i, j, k, n, rank, nullity;
    slong * pivots;
    slong * nonpivots;
    nmod_poly_mat_t tmp;
    nmod_poly_t den;

    n = mat->c;

    nmod_poly_init(den, nmod_poly_mat_modulus(mat));
    nmod_poly_mat_init_set(tmp, mat);
    rank = nmod_poly_mat_rref(tmp, den, tmp);
    nullity = n - rank;

    nmod_poly_mat_zero(res);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            nmod_poly_one(res->rows[i] + i);
    }
    else if (nullity)
    {
        pivots = flint_malloc(rank * sizeof(slong));
        nonpivots = flint_malloc(nullity * sizeof(slong));

        for (i = j = k = 0; i < rank; i++)
        {
            while (nmod_poly_is_zero(tmp->rows[i] + j))
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

        nmod_poly_set(den, tmp->rows[0] + pivots[0]);

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
                nmod_poly_set(res->rows[pivots[j]] + i,
                    tmp->rows[j] + nonpivots[i]);
            nmod_poly_neg(res->rows[nonpivots[i]] + i, den);
        }

        flint_free(pivots);
        flint_free(nonpivots);
    }

    nmod_poly_clear(den);
    nmod_poly_mat_clear(tmp);
    return nullity;
}
