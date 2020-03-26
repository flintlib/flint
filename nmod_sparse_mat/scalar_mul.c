/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

void
nmod_sparse_mat_scalar_mul(nmod_sparse_mat_t B, const nmod_sparse_mat_t A, mp_limb_t c)
{
    if (c == UWORD(0))
    {
        nmod_sparse_mat_zero(B);
    }
    else if(c == UWORD(1))
    {
        nmod_sparse_mat_set(B, A);
    }
    else if(c == B->mod.n - UWORD(1))
    {
        nmod_sparse_mat_neg(B, A);
    }
    else
    {
        slong i;
        nmod_sparse_mat_set(B, A);
        if (A->nnz > 10 && A->mod.n < UWORD_HALF)
        {
            mp_limb_t w_pr = n_mulmod_precomp_shoup(c, A->mod.n);
            for (i = 0; i<B->nnz; i++)
                B->entries[i].val = n_mulmod_shoup(c, B->entries[i].val, w_pr, A->mod.n);
        }
        else
        {
            for (i = 0; i<B->nnz; i++)
                B->entries[i].val = nmod_mul(c, B->entries[i].val, B->mod);
        }
    }
}
