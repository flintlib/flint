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

mp_limb_t nmod_sparse_vec_dot(const nmod_sparse_vec_t u, const nmod_sparse_vec_t v, nmod_t mod)
{
    slong i, j;
    mp_limb_t ret = UWORD(0);
    for (i = j = 0; i < u->nnz && j < v->nnz; )
    {
        if (u->entries[i].ind == v->entries[j].ind)
        {
            ret = nmod_addmul(ret, u->entries[i].val, v->entries[j].val, mod);
            ++i; ++j;
        }
        else if (u->entries[i].ind < v->entries[j].ind) ++i;
        else if (u->entries[i].ind > v->entries[j].ind) ++j;
    }
    return ret;
}
