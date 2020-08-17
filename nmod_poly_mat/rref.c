/*
    Copyright (C) 2011-2012 Fredrik Johansson

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
nmod_poly_mat_rref(nmod_poly_mat_t R, nmod_poly_t den, const nmod_poly_mat_t A)
{
    slong i, j, k, m, n, rank;
    slong *pivots, *nonpivots;

    rank = nmod_poly_mat_fflu(R, den, NULL, A, 0);
    m = nmod_poly_mat_nrows(R);
    n = nmod_poly_mat_ncols(R);

    /* clear bottom */
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            nmod_poly_zero(nmod_poly_mat_entry(R, i, j));

    /* Convert row echelon form to reduced row echelon form */
    if (rank > 1)
    {
        nmod_poly_t tmp, tmp2;
        nmod_poly_init(tmp, nmod_poly_mat_modulus(R));
        nmod_poly_init(tmp2, nmod_poly_mat_modulus(R));

        pivots = flint_malloc(sizeof(slong) * n);
        nonpivots = pivots + rank;

        /* find pivot positions */
        for (i = j = k = 0; i < rank; i++)
        {
            while (nmod_poly_is_zero(nmod_poly_mat_entry(R, i, j)))
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

        for (k = 0; k < n - rank; k++)
        {
            for (i = rank - 2; i >= 0; i--)
            {
                nmod_poly_mul(tmp, den, nmod_poly_mat_entry(R, i, nonpivots[k]));

                for (j = i + 1; j < rank; j++)
                {
                    nmod_poly_mul(tmp2, nmod_poly_mat_entry(R, i, pivots[j]),
                        nmod_poly_mat_entry(R, j, nonpivots[k]));
                    nmod_poly_sub(tmp, tmp, tmp2);
                }

                nmod_poly_div(nmod_poly_mat_entry(R, i, nonpivots[k]),
                    tmp, nmod_poly_mat_entry(R, i, pivots[i]));
            }
        }

        /* clear pivot columns */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    nmod_poly_set(nmod_poly_mat_entry(R, j, pivots[i]), den);
                else
                    nmod_poly_zero(nmod_poly_mat_entry(R, j, pivots[i]));
            }
        }

        flint_free(pivots);
        nmod_poly_clear(tmp);
        nmod_poly_clear(tmp2);
    }

    return rank;
}

