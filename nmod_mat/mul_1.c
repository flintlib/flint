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
_nmod_mat_mul_1(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long i, j, k;

    register mp_limb_t s;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s = 0UL;

            for (k = 0; k < A->c; k++)
            {
                s += A->rows[i][k] * B->rows[k][j];
            }

            C->rows[i][j] = s;
        }
    }

    _nmod_vec_reduce(C->entries, C->entries, A->r * B->c, C->mod);
}

void
_nmod_mat_mul_transpose_1(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long i, j, k;

    register mp_limb_t s;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->r; j++)
        {
            s = 0UL;

            for (k = 0; k < A->c; k++)
            {
                s += A->rows[i][k] * B->rows[j][k];
            }

            C->rows[i][j] = s;
        }
    }

    _nmod_vec_reduce(C->entries, C->entries, A->r * B->r, C->mod);
}
