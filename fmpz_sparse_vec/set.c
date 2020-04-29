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

void fmpz_sparse_vec_set(fmpz_sparse_vec_t dst, const fmpz_sparse_vec_t src, slong ioff) 
{
    slong i;
    if (dst == src) return;
    _fmpz_sparse_vec_resize(dst, src->nnz);
    for (i = 0; i < dst->nnz; ++i)
    {
        dst->entries[i].ind = src->entries[i].ind - ioff;
        fmpz_set(dst->entries[i].val, src->entries[i].val);
    }
}