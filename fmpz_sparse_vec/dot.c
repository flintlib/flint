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

void fmpz_sparse_vec_dot(fmpz_t ret, const fmpz_sparse_vec_t u, const fmpz_sparse_vec_t v)
{
    slong i, j;
    fmpz_zero(ret);
    for (i = j = 0; i < u->nnz && j < v->nnz; )
    {
        if (u->entries[i].ind == v->entries[j].ind)
        {
            fmpz_addmul(ret, u->entries[i].val, v->entries[j].val);
            ++i; ++j;
        }
        else if (u->entries[i].ind < v->entries[j].ind) ++i;
        else if (u->entries[i].ind > v->entries[j].ind) ++j;
    }
}
