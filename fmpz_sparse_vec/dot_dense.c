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

void fmpz_sparse_vec_dot_dense(fmpz_t ret, const fmpz_sparse_vec_t u, const fmpz *v) 
{
    slong i;
    fmpz_zero(ret);
    for (i = 0; i < u->nnz; ++i) fmpz_addmul(ret, u->entries[i].val, &v[u->entries[i].ind]);

}
