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
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"

void
nmod_mat_addmul(nmod_mat_t D, const nmod_mat_t C,
                                const nmod_mat_t A, const nmod_mat_t B)
{
    len_t m, k, n;

    m = A->r;
    k = A->c;
    n = B->c;

    if (m < NMOD_MAT_MUL_STRASSEN_CUTOFF ||
        n < NMOD_MAT_MUL_STRASSEN_CUTOFF ||
        k < NMOD_MAT_MUL_STRASSEN_CUTOFF)
    {
        _nmod_mat_mul_classical(D, C, A, B, 1);
    }
    else
    {
        nmod_mat_t tmp;
        nmod_mat_init(tmp, m, n, A->mod.n);
        nmod_mat_mul_strassen(tmp, A, B);
        nmod_mat_add(D, C, tmp);
        nmod_mat_clear(tmp);
    }
}
