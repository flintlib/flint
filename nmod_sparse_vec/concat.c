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

void nmod_sparse_vec_concat(nmod_sparse_vec_t res, const nmod_sparse_vec_t vec1,  const nmod_sparse_vec_t vec2, slong len1) 
{
    res->nnz = vec1->nnz+vec2->nnz;
    if (res->nnz == 0) nmod_sparse_vec_clear(res);
    else 
    {
        slong i;
        res->entries = flint_realloc(res->entries, res->nnz*sizeof(*res->entries));
        memcpy(res->entries, vec1->entries, vec1->nnz*sizeof(*res->entries));
        memcpy(res->entries+vec1->nnz, vec2->entries, vec2->nnz*sizeof(*res->entries));
        for (i = vec1->nnz; i < res->nnz; ++i) res->entries[i].ind += len1;
    }
}
