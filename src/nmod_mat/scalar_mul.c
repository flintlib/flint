/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
