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

void fmpz_sparse_vec_split(fmpz_sparse_vec_t res1, fmpz_sparse_vec_t res2, const fmpz_sparse_vec_t vec, slong ind)
{
    slong i, nnz1;
    fmpz_sparse_entry_struct *e1, *e2, *e;
    for (nnz1 = 0; nnz1 < vec->nnz; ++nnz1) if (vec->entries[nnz1].ind >= ind) break;

    _fmpz_sparse_vec_resize(res1, nnz1);
    _fmpz_sparse_vec_resize(res2, vec->nnz - nnz1);
    e1 = res1->entries, e2 = res2->entries;
    for (i = 0; i < vec->nnz; ++i)
    {
        e = (i < nnz1) ? e1++ : e2++;
        e->ind = vec->entries[i].ind - ((i < nnz1) ?  0 : ind);
        fmpz_set(e->val, vec->entries[i].val);
    }
}
