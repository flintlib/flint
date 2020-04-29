/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_sparse_mat.h"
#include "hashmap.h"
#include "longlong.h"

void _fmpz_sparse_mat_with_transpose_init(fmpz_sparse_mat_with_transpose_t MT, fmpz_sparse_mat_t M)
{
    slong r, c, j;
    MT->M = M;
        
    /* Construct virtual transpose */
    MT->cols = flint_calloc(M->c, sizeof(*MT->cols));
    for (r = 0; r < M->r; ++r)
        for (j = 0; j < M->rows[r].nnz; ++j)
            if (M->rows[r].entries[j].ind < M->c)
                MT->cols[M->rows[r].entries[j].ind].num++;
    for (c = 0; c < M->c; ++c)
        hashmap_init(&MT->cols[c], MT->cols[c].num);
    for (r = 0; r < M->r; ++r)
        for (j = 0; j < M->rows[r].nnz; ++j)
            if (M->rows[r].entries[j].ind < M->c)
                hashmap_put(&MT->cols[M->rows[r].entries[j].ind], r, &M->rows[r].entries[j].val);
}
