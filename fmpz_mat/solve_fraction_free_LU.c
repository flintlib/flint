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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"


void
fmpz_mat_solve_fraction_free_LU(fmpz * x, fmpz_t den, const fmpz_mat_t A,
    const fmpz * b)
{
    long dim, rank, i;
    fmpz_mat_t T;
    fmpz * tmp;
    long * perm;

    dim = A->r;

    if (dim < 1)
    {
        fmpz_set_ui(den, 1UL);
        return;
    }

    /* Compute LU decomposition in a temporary matrix */
    fmpz_mat_init_set(T, A);
    perm = malloc(dim * sizeof(long));

    rank = _fmpz_mat_rowreduce(perm, T, ROWREDUCE_FAST_ABORT);

    if (FLINT_ABS(rank) == dim)
    {
        fmpz_set(den, T->rows[dim-1] + (dim-1));
        tmp = _fmpz_vec_init(dim);
        for (i = 0; i < dim; i++)
            fmpz_set(tmp + i, b + perm[i]);

        _fmpz_mat_solve_fraction_free_LU_precomp(tmp, T);

        for (i = 0; i < dim; i++)
            fmpz_set(x + i, tmp + i);
        _fmpz_vec_clear(tmp, dim);
    }
    else
    {
        fmpz_zero(den);
    }

    free(perm);
    fmpz_mat_clear(T);
}
