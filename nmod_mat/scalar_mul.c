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
#include <gmp.h>
#include "flint.h"
#include "nmod_mat.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

void
nmod_mat_scalar_mul(nmod_mat_t B, const nmod_mat_t A, mp_limb_t c)
{
    if (c == UWORD(0))
    {
        nmod_mat_zero(B);
    }
    else if (c == UWORD(1))
    {
        nmod_mat_set(B, A);
    }
    else if (c == A->mod.n - UWORD(1))
    {
        nmod_mat_neg(B, A);
    }
    else if (A->r * A->c > 10 && A->mod.n < UWORD_HALF)
    {
        slong i, j;
        mp_limb_t w_pr = n_mulmod_precomp_shoup(c, A->mod.n);

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_mat_entry(B, i, j) = n_mulmod_shoup(
                    c, nmod_mat_entry(A, i, j), w_pr, A->mod.n);
    }
    else
    {
        slong i, j;

        for (i = 0; i < A->r; i++)
            for (j = 0; j < A->c; j++)
                nmod_mat_entry(B, i, j) = n_mulmod2_preinv(
                    nmod_mat_entry(A, i, j), c, A->mod.n, A->mod.ninv);
    }
}
