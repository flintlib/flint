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
#include "longlong.h"

void
_nmod_mat_mul_blocked_2b(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B, long block_size)
{
    long i, j, k, ii, jj, kk;

    mp_limb_t s0, s1;
    mp_limb_t t0, t1;

    nmod_mat_t S1;

    mpn_zero(C->entries, A->r * B->c);

    nmod_mat_init(S1, A->r, B->c, C->mod.n);

    for (i = 0; i < A->r; i += block_size)
    {
        for (j = 0; j < B->c; j += block_size)
        {
            for (k = 0; k < A->c; k += block_size)
            {
                for (ii = i; ii < FLINT_MIN(A->r, i + block_size); ii++)
                {
                    for (jj = j; jj < FLINT_MIN(B->c, j + block_size); jj++)
                    {
                        s0 = C->rows[ii][jj];
                        s1 = S1->rows[ii][jj];

                        for (kk = k; kk < FLINT_MIN(A->c, k + block_size); kk++)
                        {
                            umul_ppmm(t1, t0, A->rows[ii][kk], B->rows[kk][jj]);
                            add_ssaaaa(s1, s0, s1, s0, t1, t0);
                        }

                        C->rows[ii][jj] = s0;
                        S1->rows[ii][jj] = s1;
                    }
                }
            }
        }
    }

    for (i = 0; i < A->r * B->c; i++)
    {
        s0 = C->entries[i];
        s1 = S1->entries[i];

        NMOD_RED2(s0, s1, s0, C->mod);
        C->entries[i] = s0;
    }

    nmod_mat_clear(S1);
}
