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

void nmod_sparse_vec_from_entries(nmod_sparse_vec_t vec, slong * inds, mp_limb_t * vals, slong nnz) 
{
    if (nnz == 0) nmod_sparse_vec_clear(vec);
    else {
        slong i;
        vec->nnz = nnz;
        vec->entries = flint_realloc(vec->entries, nnz*sizeof(*vec->entries));
        for (i = 0; i < nnz; ++i)
        {
            vec->entries[i].ind = inds[i];
            vec->entries[i].val = vals[i];
        }
    }
}
