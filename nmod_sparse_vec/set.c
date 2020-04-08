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

void nmod_sparse_vec_set(nmod_sparse_vec_t vec, const nmod_sparse_vec_t src, slong ioff) 
{
    slong i;
    if (vec==src) return;
    if (src->nnz == 0) nmod_sparse_vec_clear(vec);
    else 
    {
        vec->entries = flint_realloc(vec->entries, src->nnz*sizeof(*vec->entries));
        memcpy(vec->entries, src->entries, src->nnz*sizeof(*vec->entries));
        vec->nnz = src->nnz;
        for (i = 0; i < vec->nnz; ++i)
        {
            vec->entries[i].ind -= ioff;
        }
    }
}