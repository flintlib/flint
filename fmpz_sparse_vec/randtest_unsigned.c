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

void fmpz_sparse_vec_randtest_unsigned(fmpz_sparse_vec_t vec, flint_rand_t state, slong nnz, slong len, flint_bitcnt_t bits)
{
    slong i, j;
    nnz = FLINT_MIN(nnz, len);
    _fmpz_sparse_vec_resize(vec, nnz);
    for (i = 0; i < nnz; ++i)
    {
        vec->entries[i].ind = i;
        do fmpz_randtest_unsigned(vec->entries[i].val, state, bits);
        while (fmpz_is_zero(vec->entries[i].val));
    }

    /* Use resevoir sampling to get random support */
    for (j = nnz; j < len; ++j)
    {
        i = n_randint(state, j+1);
        if (i < nnz) vec->entries[i].ind = j;
    }
    qsort(vec->entries, nnz, sizeof(*vec->entries), fmpz_sparse_entry_cmp);
}
