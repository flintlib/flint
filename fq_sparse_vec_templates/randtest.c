/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/
#ifdef T

#include "templates.h"

void TEMPLATE(T, sparse_vec_randtest)(TEMPLATE(T, sparse_vec_t) vec, flint_rand_t state, slong nnz, slong len, const TEMPLATE(T, ctx_t) ctx)
{
    slong i, j;
    nnz = FLINT_MIN(nnz, len);
    _TEMPLATE(T, sparse_vec_resize) (vec, nnz, ctx);
    for (i = 0; i < nnz; ++i)
    {
        vec->entries[i].ind = i;
        do TEMPLATE(T, randtest) (vec->entries[i].val, state, ctx);
        while (TEMPLATE(T, is_zero) (vec->entries[i].val, ctx));
    }

    /* Use resevoir sampling to randomize support */
    for (j = nnz; j < len; ++j)
        if ((i = n_randint(state, j+1)) < nnz) vec->entries[i].ind = j;
    qsort(vec->entries, nnz, sizeof(*vec->entries), TEMPLATE(T, sparse_entry_cmp));
}


#endif

