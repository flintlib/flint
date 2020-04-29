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

int nmod_sparse_vec_equal(const nmod_sparse_vec_t vec1, const nmod_sparse_vec_t vec2, slong ioff) 
{
    slong i;
    if (vec1->nnz != vec2->nnz) return 0;
    for (i = 0; i < vec1->nnz; ++i)
    {
        if ((vec1->entries[i].ind != vec2->entries[i].ind + ioff) || 
           (vec1->entries[i].val != vec2->entries[i].val)) return 0;
    }
    return 1;
}
