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
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long m, k, n, i, j;
    int nlimbs;

    m = A->r;
    k = A->c;
    n = B->c;

    if (k == 0)
    {
        nmod_mat_zero(C);
        return;
    }

    nlimbs = _nmod_vec_dot_bound_limbs(k, A->mod);

    if (m < NMOD_MAT_MUL_TRANSPOSE_CUTOFF ||
        n < NMOD_MAT_MUL_TRANSPOSE_CUTOFF ||
        k < NMOD_MAT_MUL_TRANSPOSE_CUTOFF)
    {
        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                nmod_mat_entry(C, i, j) = _nmod_vec_dot_ptr(A->rows[i],
                    B->rows, j, k, C->mod, nlimbs);
    }
    else
    {
        mp_ptr tmp = flint_malloc(sizeof(mp_limb_t) * k * n);

        for (i = 0; i < k; i++)
            for (j = 0; j < n; j++)
                tmp[j*k + i] = B->rows[i][j];

        for (i = 0; i < m; i++)
            for (j = 0; j < n; j++)
                nmod_mat_entry(C, i, j) = _nmod_vec_dot(A->rows[i],
                    tmp + j*k, k, C->mod, nlimbs);

        flint_free(tmp);
    }
}
