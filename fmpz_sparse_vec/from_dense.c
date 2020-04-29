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

void fmpz_sparse_vec_from_dense(fmpz_sparse_vec_t dst, const fmpz *src, slong len)
{
    slong i, nnz = 0;
    fmpz_sparse_entry_struct *e;
    for (i = 0; i < len; ++i)
        if (!fmpz_is_zero(&src[i])) ++nnz;
    _fmpz_sparse_vec_resize(dst, nnz);
    e = dst->entries;
    for (i = 0; i < len; ++i)
        if (!fmpz_is_zero(&src[i]))
            e->ind = i, fmpz_set(e->val, &src[i]), ++e;
}