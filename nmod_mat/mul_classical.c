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

/* TODO: tune for nonsquare matrices */
#define TRANSPOSE_CUTOFF 24

void
nmod_mat_mul_classical(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    mp_bitcnt_t amag;
    mp_bitcnt_t bmag;
    mp_bitcnt_t cmag;

    long m, k, n;

    int limbs;

    m = A->r;
    k = A->c;
    n = B->c;

    amag = FLINT_BIT_COUNT(A->mod.n - 1);
    bmag = FLINT_BIT_COUNT(B->mod.n - 1);
    cmag = FLINT_BIT_COUNT(k) + amag + bmag;

    limbs = cmag < FLINT_BITS ? 1 : (cmag < 2*FLINT_BITS ? 2 : 3);

    if (m < TRANSPOSE_CUTOFF || n < TRANSPOSE_CUTOFF || k < TRANSPOSE_CUTOFF)
    {
        if      (limbs == 1) _nmod_mat_mul_1(C, A, B);
        else if (limbs == 2) _nmod_mat_mul_2(C, A, B);
        else                 _nmod_mat_mul_3(C, A, B);
    }
    else
    {
        nmod_mat_t T;
        nmod_mat_init(T, B->c, B->r, B->mod.n);
        nmod_mat_transpose(T, B);

        if      (limbs == 1) _nmod_mat_mul_transpose_1(C, A, T);
        else if (limbs == 2) _nmod_mat_mul_transpose_2(C, A, T);
        else                 _nmod_mat_mul_transpose_3(C, A, T);

        nmod_mat_clear(T);
    }
}
