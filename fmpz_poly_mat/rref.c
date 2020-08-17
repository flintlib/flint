/*
    Copyright (C) 2011-2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

slong
fmpz_poly_mat_rref(fmpz_poly_mat_t R, fmpz_poly_t den, const fmpz_poly_mat_t A)
{
    slong i, j, k, m, n, rank;
    slong *pivots, *nonpivots;

    rank = fmpz_poly_mat_fflu(R, den, NULL, A, 0);
    m = fmpz_poly_mat_nrows(R);
    n = fmpz_poly_mat_ncols(R);

    /* clear bottom */
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_poly_zero(fmpz_poly_mat_entry(R, i, j));

    /* Convert row echelon form to reduced row echelon form */
    if (rank > 1)
    {
        fmpz_poly_t tmp, tmp2;
        fmpz_poly_init(tmp);
        fmpz_poly_init(tmp2);

        pivots = flint_malloc(sizeof(slong) * n);
        nonpivots = pivots + rank;

        /* find pivot positions */
        for (i = j = k = 0; i < rank; i++)
        {
            while (fmpz_poly_is_zero(fmpz_poly_mat_entry(R, i, j)))
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
                fmpz_poly_mul(tmp, den, fmpz_poly_mat_entry(R, i, nonpivots[k]));

                for (j = i + 1; j < rank; j++)
                {
                    fmpz_poly_mul(tmp2, fmpz_poly_mat_entry(R, i, pivots[j]),
                        fmpz_poly_mat_entry(R, j, nonpivots[k]));
                    fmpz_poly_sub(tmp, tmp, tmp2);
                }

                fmpz_poly_div(fmpz_poly_mat_entry(R, i, nonpivots[k]),
                    tmp, fmpz_poly_mat_entry(R, i, pivots[i]));
            }
        }

        /* clear pivot columns */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    fmpz_poly_set(fmpz_poly_mat_entry(R, j, pivots[i]), den);
                else
                    fmpz_poly_zero(fmpz_poly_mat_entry(R, j, pivots[i]));
            }
        }

        flint_free(pivots);
        fmpz_poly_clear(tmp);
        fmpz_poly_clear(tmp2);
    }

    return rank;
}

