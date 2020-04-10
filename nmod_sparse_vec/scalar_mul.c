/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"

void nmod_sparse_vec_scalar_mul_nmod(nmod_sparse_vec_t v, const nmod_sparse_vec_t u, mp_limb_t c, nmod_t mod)
{
    if (c == UWORD(0)) nmod_sparse_vec_zero(v);
    else if (c == UWORD(1)) nmod_sparse_vec_set(v, u, 0);
    else if (c==mod.n-UWORD(1)) nmod_sparse_vec_neg(v, u, mod);
    else 
    {
        slong i;
        nmod_sparse_vec_set(v, u, 0);
        for (i = 0; i < v->nnz; ++i) v->entries[i].val = nmod_mul(v->entries[i].val, c, mod);
    }
}