/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"
#include "ulong_extras.h"

void
nmod_sparse_vec_print_pretty(const nmod_sparse_vec_t vec, slong ioff, slong maxi, nmod_t mod)
{
    slong i;
    char ind_fmt[FLINT_BITS + 5];
    char val_fmt[FLINT_BITS + 5];

    flint_sprintf(ind_fmt, "%%%dwd:", n_sizeinbase(maxi, 10));
    flint_sprintf(val_fmt, "%%%dwd", n_sizeinbase(mod.n, 10));

    flint_printf("[");
    for (i = 0; i < vec->nnz; i++)
    {
        flint_printf(ind_fmt, vec->entries[i].ind - ioff);
        flint_printf(val_fmt, vec->entries[i].val);
        if (i + 1 < vec->nnz) flint_printf(" ");
    }
    flint_printf("]\n");
}

