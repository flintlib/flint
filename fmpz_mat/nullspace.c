/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_vec.h"


len_t
fmpz_mat_nullspace(fmpz_mat_t res, const fmpz_mat_t mat)
{
    len_t i, j, k, n, rank, nullity;
    len_t * pivots;
    len_t * nonpivots;
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
        pivots = flint_malloc(rank * sizeof(len_t));
        nonpivots = flint_malloc(nullity * sizeof(len_t));

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
