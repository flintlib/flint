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
#include "fmpz_poly.h"
#include "fmpz_poly_mat.h"

len_t
fmpz_poly_mat_nullspace(fmpz_poly_mat_t res, const fmpz_poly_mat_t mat)
{
    len_t i, j, k, n, rank, nullity;
    len_t * pivots;
    len_t * nonpivots;
    fmpz_poly_mat_t tmp;
    fmpz_poly_t den;

    n = mat->c;

    fmpz_poly_init(den);
    fmpz_poly_mat_init_set(tmp, mat);
    rank = fmpz_poly_mat_rref(tmp, den, tmp);
    nullity = n - rank;

    fmpz_poly_mat_zero(res);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            fmpz_poly_set_ui(res->rows[i] + i, 1UL);
    }
    else if (nullity)
    {
        pivots = flint_malloc(rank * sizeof(len_t));
        nonpivots = flint_malloc(nullity * sizeof(len_t));

        for (i = j = k = 0; i < rank; i++)
        {
            while (fmpz_poly_is_zero(tmp->rows[i] + j))
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

        fmpz_poly_set(den, tmp->rows[0] + pivots[0]);

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
                fmpz_poly_set(res->rows[pivots[j]] + i,
                    tmp->rows[j] + nonpivots[i]);
            fmpz_poly_neg(res->rows[nonpivots[i]] + i, den);
        }

        flint_free(pivots);
        flint_free(nonpivots);
    }

    fmpz_poly_clear(den);
    fmpz_poly_mat_clear(tmp);
    return nullity;
}
