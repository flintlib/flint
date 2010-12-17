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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"



void
fmpz_mat_solve_mat(fmpz_mat_t X, fmpz_t den, const fmpz_mat_t A,
                                             const fmpz_mat_t B)
{
    long i, dim, eq, equations, rank;
    fmpz_mat_t T;
    fmpz * tmp;

    FMPZ_MAT_ASSERT(A->r == A->c, "fmpz_mat_solve_mat: "
        "system matrix must be square");
    FMPZ_MAT_ASSERT(A->r == B->r, "fmpz_mat_solve_mat: "
        "different number of rows in system matrix and r.h.s.");
    FMPZ_MAT_ASSERT(A->r == X->r, "fmpz_mat_solve_mat: "
        "different number of rows in system matrix and l.h.s.");
    FMPZ_MAT_ASSERT(X->c == B->c, "fmpz_mat_solve_mat: "
        "different number of columns in l.h.s and r.h.s.");

    dim = A->r;
    equations = B->c;

    /* Potentially faster when only solving a single equation */
    if (equations <= 1 || dim < 1)
    {
        if (equations == 1)
            fmpz_mat_solve(X->entries, den, A, B->entries);
        /* else nothing to do */
        return;
    }

    /* Compute fraction-free LU decomposition */
    fmpz_mat_init_set(T, A);
    rank = _fmpz_mat_rowreduce(T, ROWREDUCE_FAST_ABORT);

    if (FLINT_ABS(rank) == dim)
    {
        long * order;
        fmpz_set(den, &T->rows[dim-1][dim-1]);

        /* Pivot order */
        order = malloc(sizeof(long) * dim);
        for (i = 0; i < dim; i++)
            order[(T->rows[i] - T->entries) / dim] = i;

        /* Solve for each column of B using precomputed LU decomposition */
        tmp = _fmpz_vec_init(dim);
        for (eq = 0; eq < equations; eq++)
        {
            for (i = 0; i < dim; i++)
                fmpz_set(&tmp[order[i]], &B->rows[i][eq]);

            _fmpz_mat_solve_fflu_precomp(tmp, T->rows, dim);

            for (i = 0; i < dim; i++)
                fmpz_set(&X->rows[i][eq], &tmp[i]);
        }
        _fmpz_vec_clear(tmp, dim);
        free(order);
    }
    else
    {
        fmpz_zero(den);
    }

    fmpz_mat_clear(T);
}
