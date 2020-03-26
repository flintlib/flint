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
#include "nmod_vec.h"

void
nmod_sparse_mat_neg(nmod_sparse_mat_t B, const nmod_sparse_mat_t A)
{
    slong i;

    nmod_sparse_mat_set(B, A);
    for (i = 0; i < A->nnz; i++)
        B->entries[i].val = nmod_neg(B->entries[i].val, B->mod);
}
