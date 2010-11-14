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

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_mat.h"


int
nmod_mat_solve_mat(nmod_mat_t X, const nmod_mat_t A, const nmod_mat_t B)
{
    long i, dim, eq, equations, rank;
    nmod_mat_t T;
    mp_limb_t * tmp;
    int result;

    dim = A->r;
    equations = B->c;

    /* Potentially faster when only solving a single equation */
    if (equations <= 1 || dim < 1)
    {
        if (equations == 1)
            return nmod_mat_solve(X->entries, A, B->entries);
        /* else nothing to do */
        return 1;
    }

    /* Compute LU decomposition */
    nmod_mat_init_set(T, A);
    rank = _nmod_mat_rowreduce(T->rows, dim, dim, ROWREDUCE_FAST_ABORT, T->mod);

    result = 0;
    if (FLINT_ABS(rank) == dim)
    {
        long * order;

        /* Pivot order */
        order = malloc(sizeof(long) * dim);
        for (i = 0; i < dim; i++)
            order[(T->rows[i] - T->entries) / dim] = i;

        /* Solve for each column of B using precomputed LU decomposition */
        tmp = nmod_vec_init(dim);
        for (eq = 0; eq < equations; eq++)
        {
            for (i = 0; i < dim; i++)
                tmp[order[i]] = B->rows[i][eq];

            _nmod_mat_solve_lu_precomp(tmp, T->rows, dim, T->mod);

            for (i = 0; i < dim; i++)
                X->rows[i][eq] = tmp[i];
        }
        nmod_vec_free(tmp);
        free(order);
        result = 1;
    }

    nmod_mat_clear(T);
    return result;
}
