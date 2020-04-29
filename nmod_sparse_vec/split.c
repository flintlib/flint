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

void nmod_sparse_vec_split(nmod_sparse_vec_t res1, nmod_sparse_vec_t res2, const nmod_sparse_vec_t vec, slong ind)
{
    slong i;
    for (i = 0; i < vec->nnz; ++i) if (vec->entries[i].ind >= ind) break;
    if (i==0) nmod_sparse_vec_clear(res1);
    else {
        res1->nnz = i;
        res1->entries = flint_realloc(res1->entries, res1->nnz*sizeof(*res1->entries));
        memcpy(res1->entries, vec->entries, res1->nnz*sizeof(*res1->entries));
    }
    if (i==vec->nnz) nmod_sparse_vec_clear(res2);
    else {
        res2->nnz = vec->nnz - i;
        res2->entries = flint_realloc(res2->entries, res2->nnz*sizeof(*res2->entries));
        memcpy(res2->entries, vec->entries+i, res2->nnz*sizeof(*res2->entries));
        for (i = 0; i < res2->nnz; ++i) res2->entries[i].ind -= ind;
    }
}
