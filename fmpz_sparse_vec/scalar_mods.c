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

void fmpz_sparse_vec_scalar_mods_fmpz(fmpz_sparse_vec_t v, const fmpz_sparse_vec_t u, const fmpz_t mod)
{
    if (fmpz_is_one(mod)) fmpz_sparse_vec_zero(v);
    else 
    {
        slong i, ind;
        fmpz_sparse_vec_set(v, u, 0);
        for (i = ind = 0; i < v->nnz; ++i) 
        {
            v->entries[ind].ind = v->entries[i].ind;
            fmpz_mods(v->entries[ind].val, v->entries[i].val, mod);
            if (!fmpz_is_zero(v->entries[ind].val)) ++ind;
        }
        _fmpz_sparse_vec_resize(v, ind);
    }
}