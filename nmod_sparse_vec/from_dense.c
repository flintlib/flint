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

void nmod_sparse_vec_from_dense(nmod_sparse_vec_t vec, mp_srcptr src, slong len)
{
    slong i;
    nmod_sparse_entry_struct *e;
    if (len == 0) nmod_sparse_vec_clear(vec);
    else {
        vec->entries = flint_realloc(vec->entries, len*sizeof(*vec->entries));
        vec->nnz = 0;
        for (i = 0; i < len; ++i)
        {
            if (src[i] == UWORD(0)) continue;
            e = &vec->entries[vec->nnz++];
            e->ind = i, e->val = src[i];
        }
        if (vec->nnz == 0) nmod_sparse_vec_clear(vec);
        else vec->entries = flint_realloc(vec->entries, vec->nnz*sizeof(*vec->entries));
    }
}