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
_nmod_mat_solve_lu(mp_limb_t * x, const nmod_mat_t A, const mp_limb_t * b)
{
    long dim, i, rank;
    nmod_mat_t T;
    mp_limb_t * tmp;
    int result;

    dim = A->r;

    nmod_mat_init_set(T, A);
    rank = _nmod_mat_rowreduce(T->rows, dim, dim, 0, T->mod);

    result = 0;
    if (FLINT_ABS(rank) == dim)
    {
        tmp = nmod_vec_init(dim);

        for (i = 0; i < dim; i++)
        {
            /* Insert in same order as the pivots */
            tmp[i] = b[(T->rows[i] - T->entries) / dim];
        }

        _nmod_mat_solve_lu_precomp(tmp, T->rows, dim, T->mod);
        mpn_copyi(x, tmp, dim);
        nmod_vec_free(tmp);
        result = 1;
    }

    nmod_mat_clear(T);
    return result;
}

int
nmod_mat_solve(mp_limb_t * x, const nmod_mat_t A, const mp_limb_t * b)
{
    return (A->r) ? _nmod_mat_solve_lu(x, A, b) : 1;
}
