/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "flint.h"
#include "fmpz_sparse_mat.h"

void
fmpz_sparse_mat_transpose(fmpz_sparse_mat_t B, const fmpz_sparse_mat_t A)
{
    slong r, c, i, j, *nnz;
    fmpz_sparse_entry_struct *Ae, *Be;
    nnz = flint_calloc(A->c, sizeof(*nnz));
    /* Get number of nnzs in each column of A (thus each row of B) */
    for (c = 0; c < A->c; ++c) 
    {
        nnz[c] = 0;
    }
    for (r = 0; r < A->r; ++r) 
    {
        for (i = 0; i < A->rows[r].nnz; ++i) 
        {
            c = A->rows[r].entries[i].ind - A->c_off;
            if (c >= A->c) break;
            nnz[c]++;
        }
    }
    /* Allocate space for nnz and reset counters */
    for (c = 0; c < A->c; ++c) 
    {
        _fmpz_sparse_vec_resize(&B->rows[c], nnz[c]);
        nnz[c] = 0;
    }
    /* Put entries into transposed matrix */
    for (r = 0; r < A->r; ++r)
    {
        for (i = 0; i < A->rows[r].nnz; ++i) 
        {
            Ae = &A->rows[r].entries[i];
            c = Ae->ind - A->c_off;
            if (c >= A->c) break;
            j = nnz[c]++;
            Be = &B->rows[c].entries[j];
            Be->ind = r;
            fmpz_set(Be->val, Ae->val);
        }
    }
    flint_free(nnz);
    B->c_off = 0;
}
