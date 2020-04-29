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

void fmpz_sparse_vec_scalar_divexact_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t c)
{
    if (fmpz_is_one(c)) fmpz_sparse_vec_set(v, u, 0);
    else if (fmpz_equal_si(c, WORD(-1))) fmpz_sparse_vec_neg(v, u);
    else 
    {
        slong i;
        fmpz_sparse_vec_set(v, u, 0);
        for (i = 0; i < v->nnz; ++i) fmpz_divexact(v->entries[i].val, v->entries[i].val, c);
    }
}