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

mp_limb_t *nmod_sparse_vec_at(nmod_sparse_vec_t vec, slong i)
{
    /* TODO: binary search */
    slong j;
    for (j = 0; j < vec->nnz; ++j)
        if (vec->entries[j].ind==i) return &vec->entries[j].val;
    return NULL;
}

