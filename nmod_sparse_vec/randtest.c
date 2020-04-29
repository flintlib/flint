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

void nmod_sparse_vec_randtest(nmod_sparse_vec_t vec, flint_rand_t state, slong nnz, slong len, nmod_t mod)
{
    slong i, j;
    mp_limb_t v;
    nnz = FLINT_MIN(nnz, len);
    vec->nnz = nnz;
    if (nnz == 0) {vec->entries = NULL; return;}

    vec->entries = flint_realloc(vec->entries, nnz*sizeof(*vec->entries));
    for (i = 0; i < nnz; ++i)
    {
        do v = n_randtest(state) % mod.n;
        while (v == UWORD(0));
        vec->entries[i].ind = i;
        vec->entries[i].val = v;
    }

    /* Use resevoir sampling to get random support */
    for (j = nnz; j < len; ++j)
    {
        i = n_randint(state, j+1);
        if (i < nnz) vec->entries[i].ind = j;
    }
    qsort(vec->entries, nnz, sizeof(*vec->entries), nmod_sparse_entry_cmp);
}
