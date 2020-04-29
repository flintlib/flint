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
#include "fmpz_sparse_vec.h"

void fmpz_sparse_vec_concat(fmpz_sparse_vec_t res, const fmpz_sparse_vec_t vec1,  const fmpz_sparse_vec_t vec2, slong len1) 
{
    slong i, nnz = vec1->nnz + vec2->nnz;
    _fmpz_sparse_vec_resize (res, nnz);
    for (i = 0; i < nnz; ++i)
    {
        fmpz_sparse_entry_struct *e = (i < vec1->nnz) ? &vec1->entries[i] : &vec2->entries[i-vec1->nnz];
        res->entries[i].ind = e->ind + ((i < vec1->nnz) ? 0 : len1);
        fmpz_set (res->entries[i].val, e->val);
    }
}
