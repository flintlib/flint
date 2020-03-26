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
#include "nmod_sparse_mat.h"
#include "ulong_extras.h"

void
nmod_sparse_mat_print_pretty(const nmod_sparse_mat_t mat)
{
    slong i, j, k=0;;
    int width;
    char row_fmt[FLINT_BITS + 5];
    char col_fmt[FLINT_BITS + 5];
    char val_fmt[FLINT_BITS + 5];

    flint_printf("<%wd x %wd sparse integer matrix mod %wu (%wd nonzeroes)>\n", 
                mat->r, mat->c, mat->mod.n, mat->nnz);

    if (!(mat->c) || !(mat->r))
        return;

    width = n_sizeinbase(mat->r, 10);
    flint_sprintf(row_fmt, "%%%dwd: [", width);

    width = n_sizeinbase(mat->c, 10);
    flint_sprintf(col_fmt, "%%%dwd:", width);

    width = n_sizeinbase(mat->mod.n, 10);
    flint_sprintf(val_fmt, "%%%dwd", width);

    for (i = 0; i < mat->r; i++)
    {
        flint_printf(row_fmt, i);

        nmod_sparse_mat_entry_struct *e = mat->entries + mat->row_starts[i];

        for (j = 0; j < mat->row_nnz[i]; j++, e++)
        {
            flint_printf(col_fmt, e->col - mat->c_off);
            flint_printf(val_fmt, e->val);
            if (j + 1 < mat->row_nnz[i])
                flint_printf(" ");
        }

        flint_printf("]\n");
    }
}

