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
_nmod_mat_mul_2(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long i, j, k;

    register mp_limb_t s0, s1;
    register mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->c; j++)
        {
            s0 = s1 = 0UL;

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, A->rows[i][k], B->rows[k][j]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }

            NMOD2_RED2(s0, s1, s0, C->mod);
            C->rows[i][j] = s0;
        }
    }
}

void
_nmod_mat_mul_transpose_2(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    long i, j, k;

    register mp_limb_t s0, s1;
    register mp_limb_t t0, t1;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < B->r; j++)
        {
            s0 = s1 = 0UL;

            for (k = 0; k < A->c; k++)
            {
                umul_ppmm(t1, t0, A->rows[i][k], B->rows[j][k]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }

            NMOD2_RED2(s0, s1, s0, C->mod);
            C->rows[i][j] = s0;
        }
    }
}
