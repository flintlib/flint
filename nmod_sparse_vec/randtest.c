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

static int nmod_sparse_entry_cmp(const void *va, const void *vb) {
    const nmod_sparse_entry_struct *a = va;
    const nmod_sparse_entry_struct *b = vb;
    if(a->ind < b->ind) return -1;
    if(b->ind < a->ind) return 1;
    return 0;
}

void nmod_sparse_vec_randtest(nmod_sparse_vec_t vec, flint_rand_t state, mp_limb_signed_t nnz, mp_limb_signed_t len, nmod_t mod) {
    nnz = FLINT_MIN(nnz, len);
    vec->nnz = nnz;
    if(nnz == 0) {vec->entries = NULL; return;}

    vec->entries = flint_realloc(vec->entries, nnz*sizeof(*vec->entries));
    slong i, j;
    mp_limb_t v;
    for(i = 0; i < nnz; ++i) {
        do v = n_randtest(state) % mod.n;
        while(v==UWORD(0));
        vec->entries[i].ind = i;
        vec->entries[i].val = v;
    }

    /* Use resevoir sampling to get random support */
    for(j = nnz; j < len; ++j) {
        i = n_randint(state, j+1);
        if(i < nnz) vec->entries[i].ind = j;
    }
    qsort(vec->entries, nnz, sizeof(*vec->entries), nmod_sparse_entry_cmp);
}
